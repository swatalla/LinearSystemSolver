# LinearSystemSolver
Exercise in developing a strongly typed linear system solver library (without using explicit loops) for my machine learning scripts.

Do not use this library - there are faster, more thorougly tested options.

Partial pivoting works, but full pivot is still under development.

## Usage:
Given a matrix of floats (64),
```fsharp
let mat10 =
    [[2; 3; 4; 5; 5; 4; 2; 7; 9; 6; 4]
     [4; 5; 2; 3; 2; 2; 6; 5; 1; 7; 6]
     [4; 2; 3; 1; 5; 5; 8; 6; 8; 1; 8]
     [1; 3; 4; 1; 2; 7; 7; 2; 3; 9; 5]
     [1; 2; 5; 2; 5; 4; 2; 1; 7; 3; 5]
     [2; 9; 7; 2; 7; 8; 6; 8; 2; 7; 7]
     [7; 7; 1; 8; 3; 5; 2; 8; 8; 4; 8]
     [1; 6; 8; 4; 7; 3; 7; 4; 3; 6; 2]
     [5; 8; 3; 7; 3; 7; 9; 6; 9; 4; 9]
     [2; 2; 7; 1; 3; 9; 9; 3; 8; 3; 3]] |> List.map (List.map float)
```
the solver may be called via
```fsharp
mat |> MatrixSystem |> solve PivotStrategy
```
where `MatrixSystem` and `PivotStrategy` are discriminated unions defined as
```fsharp
type MatrixSystem =
    | LinearSystem of float list list * float list //e.g. ax + by = cz, where ax + by and cz are separate lists
    | AugmentedMatrix of float list list //Matrix in square form, like mat10 above
    | InvalidSystem //Only used internally, will be refactored into a result type
```
and
```fsharp
type PivotStrategy =
    | PartialPivot
    | CompletePivot
```
respectively. Thus, the matrix bound to `mat10` above would be solved by calling
```fsharp
mat10 |> AugmentedMatrix |> solve PartialPivot
```
