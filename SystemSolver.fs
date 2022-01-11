namespace SystemSolver

module LinearSolver =

    type MatrixShape =
        | SquareMatrix
        | NonSquareMatrix

    type MatrixSystem =
        | LinearSystem of float list list * float list
        | AugmentedMatrix of float list list
        | InvalidSystem

    type PivotStrategy =
        | PartialPivot
        | CompletePivot

    // This should become conditioning when SVD is implemented
    type Singularity =
        | NonSingular of float list list
        | Singular

    type TriangularMatrix = 
        | LowerMatrix 
        | UpperMatrix

    type Factorization = 
        | L of float list list
        | U of float list list

    type Result<'Error, 'Solution> =
        | Error of 'Error
        | Solution of 'Solution

    // TODO: more robust pivoting
    module SystemStructure =

        // Checks that each row contains the same number of columns
        // Then, checks that the number of rows = number of columns
        let square matrix =
            matrix
            |> List.fold (fun (acc, r) row -> 
                (row |> List.fold (fun c _ -> c + 1) 0)::acc, r+1) ([], 0)
            |> fun (rows, numCols) ->
                rows
                |> List.distinct
                |> function
                    | [numRows] -> 
                        if numRows = numCols 
                            then SquareMatrix
                            else NonSquareMatrix
                    | _ -> NonSquareMatrix
            
        let augment (A: float list list) (b: float list) =
            match A |> square with
            | SquareMatrix ->
                List.map2 (fun a b -> b::(a |> List.rev) |> List.rev) A b
                |> AugmentedMatrix
            | _ -> InvalidSystem



    module Pivot =

        // Partial pivoting permutes the rows, but not the columns, such that the
        // pivot is the largest entry *in its column*

        // Complete pivoting permutes the rows and the columns such that the pivot
        // is the largest entry *in the matrix*

        // Complete pivot strategy - find row with the largest value in entire matrix,
        // bring it to the top row, first column, and permute matrix accordingly

        let inline partialPivot mat = 
            mat |> List.sortByDescending (List.head >> abs)
            //|> List.sortBy (List.head >> abs >> (~-)) // Same as sortByDescending 

        module SubMatrix =
            let inline indexPermutationByRowMaximum mat = 
                mat
                |> List.head
                |> List.rev
                |> List.tail
                |> List.rev
                |> List.indexed
                |> List.sortBy (snd >> abs)
                |> List.map fst
                |> List.toArray
                |> Array.rev

            let inline completePivotRowPermute mat =
                // sortByDescending does not seem to be stable
                let idx = mat |> indexPermutationByRowMaximum

                mat
                |> List.map (
                    List.rev
                    >> function
                        | solution::variables ->
                            solution::(variables |> List.permute (fun i -> idx.[i]))
                        | _ -> []
                    >> List.rev )

            // variable order gets messed up - need to revisit
            // need some way of keeping track of variable order
            // find some way to output the maxIdx vector
            let inline permuteBySubmatrixMaximum mat =
            
                let maxIdx =
                    mat
                    |> List.map (List.rev >> List.tail >> List.max)
                    |> List.indexed
                    |> List.sortBy (snd >> abs)
                    |> List.map fst
                    |> List.toArray
                    |> Array.rev
            
                mat
                |> List.permute (fun i -> maxIdx.[i])

            let inline reorderedVariablePermutation mat =
                mat
                |> permuteBySubmatrixMaximum
                |> indexPermutationByRowMaximum 

        let inline completePivot mat =
            match mat with
            | [_] -> mat
            | nonSingletonMatrix ->
                nonSingletonMatrix
                |> SubMatrix.permuteBySubmatrixMaximum
                |> SubMatrix.completePivotRowPermute

    module RowEchelon =

        let factor (triangle: TriangularMatrix) (subsystem: float list list) =
            subsystem
            |> List.tail
            |> List.map (fun row -> // Current row
                subsystem
                |> List.head // Pivot column
                |> fun pivot ->
                    pivot
                    |> List.map2 (fun c p ->
                        (row |> List.head)/(pivot |> List.head)
                        |> fun ratio ->
                            match triangle with
                            | LowerMatrix -> ratio
                            | UpperMatrix -> c - (ratio * p)) row)

        // partial pivot; ideally will have a full pivot that sorts columns too
        // should pivot rows, THEN pivot columns
        let reduce strategy system =
            system
            |> List.unfold (function
                | [] -> None
                | pivotStep ->
                    pivotStep
                    |> match strategy with
                        | PartialPivot -> Pivot.partialPivot 
                        | CompletePivot -> Pivot.completePivot
                    |> fun xs ->
                        (xs |> factor LowerMatrix, xs |> factor UpperMatrix)
                        |> fun (l, u) -> 
                            Some ((xs, l |> List.map List.head), u |> List.map List.tail))
            |> fun res ->
                [res |> List.map (fst >> List.head); res |> List.map snd]
                |> List.map (List.filter (List.isEmpty >> not))
                |> function
                    | a::b -> Some (a, b |> List.concat) 
                    | _ -> None

    module Conditioning =

        // Finds product of U diagonal; det(A) = det(L)det(U)... det(L) = 1
        // Should only take the U matrix
        let inline determinant u =
            u |> List.map List.head |> List.reduce (*)

        // Should only take the U matrix
        let singularity ufactor =
            ufactor
            |> determinant
            |> function
                | det when det <> 0.0 ->
                    NonSingular ufactor
                | _ -> Singular

    module Decompose =
        // only returns the upper and lower elements;
        // does not return full-rank matrices
        let lu strategy system =
            system
            |> RowEchelon.reduce strategy
            |> function
                | Some xs -> 
                    Solution (xs |> fst |> U, xs |> snd |> L)
                | None -> 
                    Error "Error reducing matrix to Row Echelon Form"
        // need a separate formatter to return full matrix, but LU decomp is
        // not main purpose of this lib
    
    module Solve =

        let backsub = function
            | [] -> None
            | system ->
                system
                |> List.head
                |> function
                    | b::Ax -> 
                        b / (Ax |> List.head)
                        |> fun soln -> 
                            system
                            |> List.map (fun row ->
                                row
                                |> List.tail
                                |> fun rs -> (row |> List.head) - (rs |> List.head) * soln::(rs |> List.tail))
                            |> fun sys -> Some (soln, sys |> List.tail)
                    | _ -> None

        // This takes the upper triangular matrix from the LU decomposition
        // Terminates with error if matrix is singular
        let rref = function
            | U upper ->
                upper
                |> Conditioning.singularity
                |> function
                    | NonSingular umatrix ->
                        umatrix
                        |> List.rev
                        |> List.map List.rev
                        |> backsub
                        |> List.unfold (function 
                            | Some eq -> Some (fst eq, eq |> snd |> backsub) 
                            | None -> None)
                        |> List.rev
                        |> Solution
                    | Singular -> Error "Singular Matrix: No unique solutions"
            | _ -> Error "RREF requires an Upper triangular matrix"

    // Gaussian Elimination via back-substitution
    let rec solve strategy = function
        | AugmentedMatrix augmentedLinearSystem ->
            augmentedLinearSystem
            |> Decompose.lu strategy
            |> function
                | Solution decomposition ->
                    decomposition
                    |> fst
                    |> Solve.rref
                    |> fun rowReducedSystem ->
                        match strategy with
                        | CompletePivot ->
                            match rowReducedSystem with
                            | Solution solution ->
                                let permutationVector = 
                                    augmentedLinearSystem 
                                    |> Pivot.SubMatrix.reorderedVariablePermutation
                            
                                solution
                                |> List.permute (fun i -> permutationVector.[i])
                                |> Solution
                            | _ -> Error "No known solutions"
                        | PartialPivot -> rowReducedSystem                    
                | Error err -> Error err
        | LinearSystem (A, b) -> // Augments the linear system and recurses
            SystemStructure.augment A b
            |> solve strategy
        | InvalidSystem ->
            Error "Solver encountered an invalid system"