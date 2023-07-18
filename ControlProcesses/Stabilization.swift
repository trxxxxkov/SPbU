/**STABILIZATION**

This program is supposed to be used in an online compiler such as onlineGDB.

BRIEF GUIDE:
 The program solves stabilization problem for 3 types of time-invariant systems,
 providing detailed information about intermediate calculations:
  1-st type: FULLY controllable systems with SCALAR input;
  2-nd type: FULLY controllable systems with VECTOR input;
  3-rd type: PARTICIALLY controllable systems with SCALAR input.
 
 
INITIALIZATION:
 Matlab-like input style is used here, which means that elements in a single
 row may or may not be separated by comma:
    
    (#1) 1 2; 3 4; 5 6 [Enter] -->  [ 1, 2;
                                      3, 4;
    (#2) 1,2; 3,4; 5,6 [Enter] -->    5, 6 ]
 
 and after each row there should be a semicolon indicating its ending(except the last one,
 where the end sign is optional):
 
    (#3) 1 2 3 4; 5 6 7 8 [Enter]  -->  [ 1, 2, 3, 4;
    (#4) 1 2 3 4; 5 6 7 8; [Enter] -->    5, 6, 7, 8 ].
 
 Square brackets are also available, but why would you use them, huh?
 (#1) is recommended as the most succinct one.
 
 Eigenvalues are entered the same way as a regular row-vector, all the aforementioned rules
 are appliable.
 
 
 KNOWN ERRORS:
 
 1) If you type any characters like a cyrillic letter, " ` " or " ' " while entering a matrix,
 you will get an error even if the character will be deleted before you press Enter.
 If this happens, just restart the program.
 
 2) Each solution ends with a check of the answer with a comparison of two
 polynomials, which are not going to be the same if any fractional numbers were involved
 in calculations. There will be a warning that answer may not be correct because of
 a rounding errors, that make polynomials *look* different. You may ignore that warning.
 
 
 OUTPUT:
 The result of the program is a step-by-step calculations for finding a matrix M,
 that allows you to find a control:
 
    u = Mx.
 
 The assistive messages in the terminal window guide through the whole process,
 though certain awareness about the methods is required. To achieve it you can use
 this manual:
 "Оптимальная стабилизация линейных систем:
 Учеб. пособие / Тамасян Г. Ш., Фоминых А. В. - СПб.:
 Изд-во ВВМ, 2022. - 66 с. Библ. 14 назв."
 ISBN 978-5-9651-1405-4
 
*/

struct Matrix {                                                                 // Efficient determinant; somewhat efficient inversion... A convenient matrix structure for daily usage
    var body: [Double] = []
    var dimension: (y: Int, x: Int)
    
    var rows: [[Double]] {                                                      // Returns the whole matrix in the format of a 2-dimensional array
        var answer: [[Double]] = Array(repeating: Array(repeating: 0, count: dimension.x), count: dimension.y)
        for i in 0..<dimension.y {
            for j in 0..<dimension.x {
                answer[i][j] = self[i, j]
            }
        }
        return answer
    }
    
    var columns: [[Double]] {                                                   // Returns the whole matrix in the format of a 2-dimensional array
        var answer: [[Double]] = Array(repeating: Array(repeating: 0, count: dimension.y), count: dimension.x)
        for j in 0..<dimension.x {
            for i in 0..<dimension.y {
                answer[j][i] = self[i, j]
            }
        }
        return answer
    }

    var characteristicPolynomial: [Double] {                                    // Uses Danilevsky method to obtain coefficients of the characteristic polynomial
        if self.body.count == 1 {
            return [1, -self[0, 0]]
        }
        
        var B = Matrix(dimension: self.dimension, body: nil)
        var D = Matrix(dimension: self.dimension, body: self.body)              // A matrix to be transformed into Frobenius normal form
        var answer: [Double] = []
        let n = self.dimension.y - 1                                            //k = (n - step + 1) - a transition function between what I'd read and what I wrote
        for step in 1...(self.dimension.y - 1) {
            
            if D[n - step + 1, n - step] == 0 {                                 // The element is going to be used as a divider. This is a try to exchange it
                for m in 0..<(n - step) {                                       // for another one at the left(with smaller second index) using determinant rules
                    if D[n - step + 1, m] != 0 {
                        D.swapColumns(n - step, m)
                        D.swapRows(n - step, m)
                        break
                    }
                }
                
                if D[n - step + 1, n - step] == 0 {                             // If there is no available elements for an exchange, an original
                                                                                // polynomial can be calculated as a product of two:
                    let det1: [Double] = D.submatrix(from: (0, 0),
                                                        to: (n - step, n - step)).characteristicPolynomial // Upper left...
                    var det2: [Double] = [1]                                                               // ...and lower right
                    for k in (n - step + 1)..<D.dimension.x {
                        det2.append(-D[(n - step + 1), k])
                    }
                    
                    var det1det2: [Double] = Array(repeating: 0, count: det1.count + det2.count - 1)        // Product of the two polynomials
                    for multiplier1 in 0..<det1.count {
                        for multiplier2 in 0..<det2.count {
                            det1det2[det1det2.count - 1 - (det1.count - 1 - multiplier1 + det2.count - 1 - multiplier2)] += det1[multiplier1] * det2[multiplier2]
                        }
                    }
                    answer = det1det2
                    return answer
                }
            }
            
            B = Matrix(dimension: self.dimension, body: nil)
            for i in 0..<B.dimension.y {                                        // A normal process with non-zero divider on the way of the Danilevsky method
                for j in 0..<B.dimension.x {
                    if i == j {
                        B[i, j] = 1
                    }
                    if i == n - step {
                        if j != n - step {
                            B[n - step, j] = -D[n - step + 1, j] / D[n - step + 1, n - step]
                        } else {
                            B[n - step, j] = 1 / D[n - step + 1, n - step]
                        }
                    }
                }
            }
            D = (B.inversed() * D) * B
        }
        answer = [1] + (-1 * D).rows[0]
        return answer
    }
    
    var rank: Int {                                                             // Calculates number of non-zero rows after forward elimination of the Gauss method
        if self.body.count == 1 {
            return 1
        }
    
                                                                                // NOTE: the problem in numerical calculation of a rank is
                                                                                // rounding errors that raise question of what is considered zero..
    
        let thisIsConsideredZero: Double = 0.00001                              // Here is my opinion: zero is anything less than 10^-5
        
        var temporaryMatrix = self
        for i in 0..<temporaryMatrix.dimension.x {                              // Zero columns can interrupt forward elimination, so they should be removed:
            if temporaryMatrix.columns[i].min() == 0 && temporaryMatrix.columns[i].max() == 0 { // Check if the column contains zeros;
                temporaryMatrix.swapColumns(i, temporaryMatrix.dimension.x - 1) // Move zero-column to the edge of the matrix;
                temporaryMatrix = temporaryMatrix.submatrix(from: (0, 0),       // Use submatrix which does not include zero-column;
                                                            to: (temporaryMatrix.dimension.y - 1, temporaryMatrix.dimension.x - 2))
                return temporaryMatrix.rank                                     // If there were several zero-columns, another check is required,
            }                                                                   // so it is easier to restart the process
        }
        
        
        var answer: Int
        if self.dimension.x < self.dimension.y {                                // The value of the answer will only be rewritten if rank < self.dimension,
            answer = self.dimension.x                                           // so it equals to a dimension of a squared matrix self
        } else {
            answer = self.dimension.y
        }
        
        for d in 0..<answer {
            if temporaryMatrix[d, d] == 0 {                                     // Exclusion zero elements from the diagonal is needed
                for i in d..<temporaryMatrix.dimension.y {
                    var rowChecker: Int = 0
                    for horizontalElementIndex in d..<temporaryMatrix.dimension.x {
                        if temporaryMatrix[i, horizontalElementIndex] != 0 {
                            rowChecker += 1
                        }
                    }
                    if rowChecker == 0 && i == (temporaryMatrix.dimension.y - 1) {   // If there is a column that has only zeros below some diagonal element
                        answer = d                                                              // and this element is also zero, rank has been found
                        return answer
                    }
                    if temporaryMatrix[i, d] != 0 {
                        temporaryMatrix.swapRows(d, i)
                        break
                    }
                }
            }
                                                                                // Forward elimination of the Gauss method
            for i in (d+1)..<temporaryMatrix.dimension.y {
                for j in stride(from: temporaryMatrix.dimension.x - 1, through: d, by: -1) {    // Element at the left participate in the formula for the other elements
                    temporaryMatrix[i, j] -= temporaryMatrix[d, j] * (temporaryMatrix[i, d] / temporaryMatrix[d, d])     // of the row so it should be changed last
                    if abs(temporaryMatrix[i, j]) < thisIsConsideredZero {
                        temporaryMatrix[i, j] = 0
                    }
                }
            }
        }
        return answer
    }
    
    init(dimension: (y: Int, x: Int), body: [Double]?) {
        if let body = body {
            self.body = body
        } else {
            self.body = Array(repeating: 0, count: dimension.x * dimension.y)
        }
        self.dimension = dimension
    }
    
    init(sBody: String) {                                                       // Removes [], if there are any, split the string by ';' into rows
        var bodyWithoutSquareBraces = sBody                                     // and split rows into elements by ' ', ',' and ', '
        bodyWithoutSquareBraces.removeAll(where: { $0 == "[" || $0 == "]" })
        var rows: [String] = []
        for row in bodyWithoutSquareBraces.split(whereSeparator: {$0 == ";"}) {
            rows.append(String(row))
            if row.firstIndex(of: ",") != nil {
                var rowWithoutSpaces = row
                rowWithoutSpaces.removeAll(where: {$0 == " "})
                for element in rowWithoutSpaces.split(whereSeparator: {$0 == ","}) {
                    body.append(Double(element)!)
                }
                
            } else {
                for element in row.split(whereSeparator: {$0 == " "}) {
                    body.append(Double(element)!)
                }
            }
        }
        dimension = (rows.count, body.count / rows.count)
    }
    
    subscript(_ y: Int, _ x: Int) -> Double {                                   // Output format: [verticalCoordinate, horizontalCoordinate]
        get {
            return body[y * dimension.x + x]
        }
        set {
            body[y * dimension.x + x] = newValue
        }
    }
    
    static func + (mLeft: Matrix, mRight: Matrix) -> Matrix {                   // Matrix + Matrix
        var answer = Matrix(dimension: mLeft.dimension, body: nil)
        for i in 0..<answer.body.count {
            answer.body[i] = mLeft.body[i] + mRight.body[i]
        }
        return answer
    }
    
    static func - (mLeft: Matrix, mRight: Matrix) -> Matrix {                   // Matrix - Matrix
        var answer = Matrix(dimension: mLeft.dimension, body: nil)
        for i in 0..<answer.body.count {
            answer.body[i] = mLeft.body[i] - mRight.body[i]
        }
        return answer
    }
    
    static func * (dLeft: Double, mRight: Matrix) -> Matrix {                   // Double * Matrix
        var answer = Matrix(dimension: mRight.dimension, body: nil)
        for i in 0..<mRight.body.count {
            answer.body[i] = dLeft * mRight.body[i]
        }
        return answer
    }
    
    static func * (mLeft: Matrix, mRight: Matrix) -> Matrix {                   // Matrix * Matrix
    
        if mLeft.body.count == 1 {                                              // If one of the operands is just a number(a single-element matrix) ->
            return mLeft[0, 0] * mRight                                         // delegate the call to an appropriate function
        }
        if mRight.body.count == 1 {
            return mRight[0, 0] * mLeft
        }
        
        var answer = Matrix(dimension: (mLeft.dimension.y, mRight.dimension.x), body: nil)
        for i in 0..<answer.dimension.y {
            for j in 0..<answer.dimension.x {
                for k in 0..<mLeft.dimension.x {
                    answer[i, j] += mLeft[i, k] * mRight[k, j]                  // Real matrix product
                }
            }
        }
        return answer
    }
    
    mutating func setRow(at rowIndex: Int, _ newRow: [Double]) {
        for j in 0..<dimension.x {
            self[rowIndex, j] = newRow[j]
        }
    }
    
    mutating func setColumn(at columnIndex: Int, _ newColumn: [Double]) {
        for i in 0..<dimension.y {
            self[i, columnIndex] = newColumn[i]
        }
    }
    
    mutating func swapRows(_ firstRowIndex: Int, _ secondRowIndex: Int) {
        let temporaryRow: [Double] = rows[firstRowIndex]
        setRow(at: firstRowIndex, rows[secondRowIndex])
        setRow(at: secondRowIndex, temporaryRow)
    }
    
    mutating func swapColumns(_ firstColumnIndex: Int, _ secondColumnIndex: Int) {
        let temporaryColumn: [Double] = columns[firstColumnIndex]
        setColumn(at: firstColumnIndex, columns[secondColumnIndex])
        setColumn(at: secondColumnIndex, temporaryColumn)
    }
    
    func submatrix(from upperLeft: (y: Int, x: Int), to lowerRight: (y: Int, x: Int)) -> Matrix {   // Returns submatrix which borders in original matrix
        var answer = Matrix(dimension: (lowerRight.y - upperLeft.y + 1, lowerRight.x - upperLeft.x + 1), body: nil) // are defined by two elements: upper left and lower right
        for i in upperLeft.y...lowerRight.y {
            for j in upperLeft.x...lowerRight.x {
                answer[i - upperLeft.y, j - upperLeft.x] = self[i, j]
            }
        }
        return answer
    }
    
    func det() -> Double {                                                      // For 1x1, 2x2 and 3x3 is calculated by explicit formulas.
        if dimension.x == 1 && dimension.y == 1 {                               // 4x4+ is calculated with modification of the Gauss method
            return self[0, 0]
        } else if dimension.x == 2 && dimension.y == 2 {
            return self[0, 0] * self[1, 1] - self[0, 1] * self[1, 0]
        } else if dimension.x == 3 && dimension.y == 3 {
            return self[0, 0] * self[1, 1] * self[2, 2] + self[1, 0] * self[2, 1] * self[0, 2] + self[2, 0] * self[0, 1] * self[1, 2] - self[2, 0] * self[1, 1] * self[0, 2] - self[1, 0] * self[0, 1] * self[2, 2] - self[0, 0] * self[2, 1] * self[1, 2]
        } else {
            var temporaryMatrix = self
            var answer: Double = 1
            for d in 0..<temporaryMatrix.dimension.x {
                if temporaryMatrix[d, d] == 0 {                                 // Check if there is a zero element to be exclude from the diagonal
                    for i in d..<temporaryMatrix.dimension.y {
                        if temporaryMatrix[i, d] == 0 && i == (temporaryMatrix.dimension.y - 1) {
                            return 0
                        }
                        if temporaryMatrix[i, d] != 0 {
                            temporaryMatrix.swapRows(d, i)
                            answer *= -1
                            break
                        }
                    }
                }
                for i in (d+1)..<temporaryMatrix.dimension.y {                  // Forward elimination from the Gauss method
                    for j in stride(from: temporaryMatrix.dimension.x - 1, through: d, by: -1) {
                        temporaryMatrix[i, j] -= temporaryMatrix[d, j] * (temporaryMatrix[i, d] / temporaryMatrix[d, d])
                    }
                }
                answer *= temporaryMatrix[d, d]
            }
            return answer
        }
    }
    
    func transposed() -> Matrix {
        var answer = Matrix(dimension: (dimension.x, dimension.y), body: nil)
        for i in 0..<dimension.y {
            for j in 0..<dimension.x {
                answer[j, i] = self[i, j]
            }
        }
        return answer
    }
    
    func inversed() -> Matrix {                                                 // 1x1, 2x2 and 3x3 are computed by explicit formulas.
        var answer = Matrix(dimension: dimension, body: nil)                    // 4x4+ inversion is based on the Cramer's rule
        if dimension.x == 1 && dimension.y == 1 {
            answer[0, 0] = 1 / self.det()
        } else if dimension.x == 2 && dimension.y == 2 {
            answer[0, 0] = self[1, 1] / self.det()
            answer[0, 1] = -self[0, 1] / self.det()
            answer[1, 0] = -self[1, 0] / self.det()
            answer[1, 1] = self[0, 0] / self.det()
        } else if dimension.x == 3 && dimension.y == 3 {
            answer[0, 0] = (self[1, 1] * self[2, 2] - self[2, 1] * self[1, 2]) / self.det()
            answer[0, 1] = (self[0, 2] * self[2, 1] - self[0, 1] * self[2, 2]) / self.det()
            answer[0, 2] = (self[0, 1] * self[1, 2] - self[0, 2] * self[1, 1]) / self.det()
            answer[1, 0] = (self[1, 2] * self[2, 0] - self[1, 0] * self[2, 2]) / self.det()
            answer[1, 1] = (self[0, 0] * self[2, 2] - self[0, 2] * self[2, 0]) / self.det()
            answer[1, 2] = (self[1, 0] * self[0, 2] - self[0, 0] * self[1, 2]) / self.det()
            answer[2, 0] = (self[1, 0] * self[2, 1] - self[2, 0] * self[1, 1]) / self.det()
            answer[2, 1] = (self[2, 0] * self[0, 1] - self[0, 0] * self[2, 1]) / self.det()
            answer[2, 2] = (self[0, 0] * self[1, 1] - self[1, 0] * self[0, 1]) / self.det()
        } else {
            for i in 0..<dimension.y {
                for j in 0..<dimension.x {
                    var submatrixBody: [Double] = []                            // A minor of a smaller dimension
                    for m in 0..<dimension.y {
                        for n in 0..<dimension.x {
                            if m != i && n != j {
                                submatrixBody.append(self[m, n])
                            }
                        }
                    }
                    answer[i, j] = (1/self.det()) * Matrix(dimension: (dimension.y - 1, dimension.x - 1), body: submatrixBody).det() * ((i + j) % 2 == 0 ? 1 : -1)
                }
            }
            answer = answer.transposed()
        }
        return answer
    }
    
    func isStable() -> Bool {                                                   // Uses Routh-Hurwitz stability criterion and characteristicPolynomial propety of the matrix structure
        var answer: Bool = true
        let polynomial: [Double] = self.characteristicPolynomial
        var HurwitzMatrix = Matrix(dimension: (polynomial.count - 1, polynomial.count - 1), body: nil)
        for i in 0..<HurwitzMatrix.dimension.y {
            for j in 0..<HurwitzMatrix.dimension.x {
                if (j - i + 1 + j) >= 0 && (j - i + 1 + j) < polynomial.count {
                    HurwitzMatrix[i, j] = polynomial[j - i + 1 + j]
                } else {
                    HurwitzMatrix[i, j] = 0
                }
            }
        }
        for d in 0..<HurwitzMatrix.dimension.y {
            if HurwitzMatrix.submatrix(from: (0, 0), to: (d, d)).det() <= 0 {
                answer = false
                break
            }
        }
        return answer
    }
    
    func printBody(withPrecision: Int?) {                                       // Note: rounding is not mathematical: it just drops all signs after the specified one
        for i in 0..<dimension.y {
            for j in 0..<dimension.x {
                var element = String(self[i, j])
                if var precision = withPrecision {
                    if precision == 0 {
                        precision = -1                                          // This erases a dot between an integer and a fractional parts of the number
                    }
                    
                    if element.contains(where: {$0 == "e" || $0 == "E"}) {      // Processing numbers in scientific format
                        var temporaryString = element
                        
                        temporaryString.remove(at: temporaryString.firstIndex(of: ".")!)
                        temporaryString = String(temporaryString[..<temporaryString.firstIndex(of: "e")!])  // Extract non-zero digits from the number
                        
                                                                                // Combine the string from the previous step with the proper amount of
                                                                                // zeros, acquired from an exponent of the original number
                        element = "0." + String(repeating: "0", count: Int(element[element.index(element.firstIndex(of: "e")!, offsetBy: 2)...])! - 1) + temporaryString
                    }
                    
                    let indexOfTheDot = element.firstIndex(of: ".")                                                                 // Precision adjustments:
                    if precision - element[element.firstIndex(of: ".")!...].count >= 0 {                                            // If string representation of the elements is shorter than
                        element += String(repeating: "0", count: precision - element[element.firstIndex(of: ".")!...].count + 1)    // required precision, it is filled with zeros
                    }
                    element = String(element[...element.index(indexOfTheDot!, offsetBy: precision)])
                }
                if i == 0 && j == 0 {
                    print("[", separator: "", terminator: "")
                } else if i != 0 && j == 0 {
                    print(" ", separator: "", terminator: "")
                }
                if j == (dimension.x - 1) {
                    print(element, separator: "", terminator: "")
                } else {
                    print(element, ",", separator: "", terminator: "")
                }
            }
            if i == (dimension.y - 1) {
                print("]\n\n")
            } else {
                print(";\n")
            }
        }
    }
    
    
}

struct Butler {                                                                 // Auxiliary structure that takes care about human-readable initialization process
                                                                                // and some other visual amenities
    let P: Matrix
    let Q: Matrix
    let targetEigenvalues: [Double]
    var S: Matrix
    
    init() {
        var stringInput: String
        
        print("P: ", terminator: "")
        stringInput = readLine()!
        P = Matrix(sBody: stringInput)
        
        print("Q: ", terminator: "")
        stringInput = readLine()!
        Q = Matrix(sBody: stringInput)
        
        S = controllabilityMatrix(P, Q)                                         // The matrix is needed here to aquire information about target eigenvalues' number
        
        print("Эталонные собственные числа(\(S.rank) шт.): ", terminator: "")
        stringInput = readLine()!
        targetEigenvalues = Matrix(sBody: stringInput).rows[0]
        
        
        if isFullyControllable(S) {
            if hasScalarInput(Q) {
                fullyControllableWithScalarInput(P, Q, S, targetEigenvalues)
            } else {
                fullyControllableWithVectorInput(P, Q, S, targetEigenvalues)
            }
        } else {
            if hasScalarInput(Q) {
                partiallyControllableWithScalarInput(P, Q, targetEigenvalues)
            } else {
                print("Алгоритм стабилизации НЕПОЛНОСТЬЮ управляемых систем с ВЕКТОРНЫМ входом в этой программе не реализован. Да и не дают у нас такие упражнения... В любом случае: желаю удачи!")
            }
        }
    }
}

func isFullyControllable(_ S: Matrix) -> Bool {
    var answer: Bool
    if S.rank == S.dimension.y {
        answer = true
    } else {
        answer = false
    }
    return answer
}

func hasScalarInput(_ Q: Matrix) -> Bool {
    var answer: Bool
    if Q.dimension.x == 1 {
        answer = true
    } else {
        answer = false
    }
    return answer
}

func controllabilityMatrix(_ P: Matrix, _ Q: Matrix) -> Matrix {                // Builds S, which is needed for Kalman rank condition of controllability for time-invariant systems
    var answer: Matrix
    if Q.dimension.x == 1 {                                                     // If control is scalar, simply uses columns Q, PQ,..., P^(n-1)Q to aquire the matrix
        answer = Matrix(dimension: P.dimension, body: nil)
        var PQ = Q
        for i in 0..<P.dimension.x {
            answer.setColumn(at: i, PQ.columns[0])
            PQ = P * PQ
        }
    } else {                                                                    // If control is a vector, all columns of Q, PQ,... are added together
        answer = Matrix(dimension: (P.dimension.y, P.dimension.x * Q.dimension.x), body: nil)
        var PQ = Q
    outerLoop: for i in stride(from: 0, to: answer.dimension.x, by: Q.dimension.x) {    // Outer loop is for setting columns in S;
        innerLoop: for QColumnNumber in 0..<Q.dimension.x {                                    // Inner loop for selecting an appropriate column of P^(i)Q to use;
                answer.setColumn(at: i + QColumnNumber, PQ.columns[QColumnNumber])
                if answer.rank == P.dimension.y {                               // If maximum rank for S has been scored, it is cutted
                    answer = answer.submatrix(from: (0, 0),                     // to a square matrix and returned immediately
                                                to: (answer.dimension.y - 1, i + QColumnNumber))
                    break outerLoop
                }
            }
            PQ = P * PQ
        }
    }
    return answer
}

func polynomialsProduct(_ p1: [Double], _ p2: [Double]) -> [Double] {           // Note: answer[0] is a coefficient of the variaable of the highest degree:
    var answer: [Double] = Array(repeating: 0, count: p1.count + p2.count - 1)  // Example: 8x^4 - 5x^3 - 113x^1 + 334 --> [8, -5, 0, -113, 334]
    for multiplier1 in 0..<p1.count {
        for multiplier2 in 0..<p2.count {
            answer[answer.count - 1 - (p1.count - 1 - multiplier1 + p2.count - 1 - multiplier2)] += p1[multiplier1] * p2[multiplier2]
        }
    }
    return answer
}

func printPolynomial(_ polynomial: [Double], rightSideOfEquation: [Double]?, singleVariable: Bool) {    // Is used for writing accompanying information about calculations into terminal window
    if singleVariable {
        for n in 0..<polynomial.count {
            if n == (polynomial.count - 1) {
                print("\(polynomial[n] >= 0 ? " +" : " ")", "\(polynomial[n])", separator: "", terminator: "")
            } else if n == 0 {
                print("\(polynomial[n])x^\(polynomial.count - 1 - n)", separator: "", terminator: "")
            } else {
                print("\(polynomial[n] >= 0 ? " +" : " ")", "\(polynomial[n])x^\(polynomial.count - 1 - n)", separator: "", terminator: "")
            }
        }
        if let rightSide = rightSideOfEquation {
            print(" = \(rightSide[0])\n")
        } else {
            print(" = 0\n")
        }
    } else {
        for n in 0..<polynomial.count {
            if n == 0 {
                print("\(polynomial[n])x\(n + 1)", separator: "", terminator: "")
            } else {
                print("\(polynomial[n] >= 0 ? " +" : " ")\(polynomial[n])x\(n + 1)", separator: "", terminator: "")
            }
        }
        if let rightSide = rightSideOfEquation {
            print(" = \(rightSide[0])\n")
        } else {
            print(" = 0\n")
        }
    }
}

func fullyControllableWithScalarInput(_ P: Matrix, _ Q: Matrix, _ S: Matrix, _ targetEigenvalues: [Double]) {
    print("Решается задача стабилизации системы со СКАЛЯРНЫМ входом для случая ПОЛНОЙ управляемости.\n")
    
    print("Построим матрицу S:")
    S.printBody(withPrecision: nil)
    print("Ранг S = n = \(P.dimension.x). Система полностью управляема.")
    
    var M: Matrix
    let polynomialOfP = P.characteristicPolynomial
    print("Найдем коэффициенты характеристического полинома матрицы P:")
    printPolynomial(polynomialOfP, rightSideOfEquation: nil, singleVariable: true)
    
    print("Строим верхнетреугольную матрицу K:")
    var K = Matrix(dimension: P.dimension, body: nil)
    for i in 0..<K.dimension.y {
        for j in 0..<K.dimension.x {
            if i == j {
                K[i, j] = 1
            } else if i < j {
                K[i, j] = polynomialOfP[j - i]
            }
        }
    }
    K.printBody(withPrecision: nil)
    
    print("Выпишем эталонный многочлен:")
    var targetPolynomial: [Double] = [1, -targetEigenvalues[0]]
    for i in 1..<targetEigenvalues.count {
        targetPolynomial = polynomialsProduct(targetPolynomial, [1, -targetEigenvalues[i]])
    }
    printPolynomial(targetPolynomial, rightSideOfEquation: nil, singleVariable: true)
    
    print("Вычислим вектор строку разностей коеффициентов хар. полиномов Г:")
    var G = Matrix(dimension: (1, P.dimension.x), body: nil)
    for i in 0...(P.dimension.x - 1) {
        G[0, i] = polynomialOfP[i + 1] - targetPolynomial[i + 1]
    }
    G.printBody(withPrecision: nil)
    
    print("K^(-1):")
    K.inversed().printBody(withPrecision: nil)
    print("S^(-1):")
    S.inversed().printBody(withPrecision: nil)
    print("Г * K^(-1):")
    (G * K.inversed()).printBody(withPrecision: nil)
    print("M = Г * K^(-1) * S^(-1):")
    M = G * K.inversed() * S.inversed()
    M.printBody(withPrecision: nil)
    
    print("Проверим найденное стабилизирующее управление. Вычислим P + (Q * M):")
    (P + Q * M).printBody(withPrecision: nil)
    print("Найдем характеристический полином P + (Q * M):")
    printPolynomial((P + Q * M).characteristicPolynomial, rightSideOfEquation: nil, singleVariable: true)
    print("Сравним его с эталонным полиномом:")
    printPolynomial(targetPolynomial, rightSideOfEquation: nil, singleVariable: true)
    if (P + Q * M).characteristicPolynomial == targetPolynomial {
        print("Отлично! Мне кажется, их корни совпадают!")
    } else {
        print("Отлично! Мне кажется, их корни не совпадают! Проверьте, не является ли это следствием погрешности округления вычислений и начните с начала.")
    }
}

func fullyControllableWithVectorInput(_ P: Matrix, _ Q: Matrix, _ S: Matrix, _ targetEigenvalues: [Double]) {
    print("Решается задача стабилизации системы с ВЕКТОРНЫМ входом для случая ПОЛНОЙ управляемости.\n")
    
    print("Построим матрицу S:")
    S.printBody(withPrecision: nil)
    print("Ранг S = n = \(P.dimension.x). Система полностью управляема.\n")
    
    var M: Matrix
    
    print("Теперь вычислим матрицу S_waved.")
    var S_waved = Matrix(dimension: P.dimension, body: nil)                     // There used another algorithm of building controllability matrix:
                                                                                // Add i'th column of Q, PQ,... then go to the next i until you get
                                                                                // a fully controllable matrix with rank of S_waved = n.
    
    var linearIndependencyCheck = 0                                             // Indicator of S_waved's rank changes
    var numberOfQColumnToUse = 0
    var whereQColumnChanged = 0
    var PQ = Q
    var S_wavedInizializationIndex: Int = 0                                     // Regular counter variable like 'i', that is used in while loop to set columns of S_waved
    
    var linearIndependentColumnsOfPQ: [Int] = []                                // Collection of numbers that indicate how much linear independent column were
                                                                                // aquired from a particular column of Q being multiplied by P, P^2, P^3,..
                                                                                // - is involved in building subsystems of P_waved matrix
    
    while S_wavedInizializationIndex < S_waved.dimension.x {                                // S_waved a priori has the rank equals to n, so upper bound for the initialization
        S_waved.setColumn(at: S_wavedInizializationIndex, PQ.columns[numberOfQColumnToUse]) // index is alredy known
        PQ = P * PQ
        if S_waved.rank == P.dimension.x {                                                              // If required amount of linear independent columns has been collected, the process
            linearIndependentColumnsOfPQ.append(S_wavedInizializationIndex - whereQColumnChanged + 1)   // is interrupted
            break
        }
        if S_waved.rank <= linearIndependencyCheck {                            // If last added column hasn't increased the rank of S_waved,
            numberOfQColumnToUse += 1                                           // next column of Q is now in use;
            
            PQ = Q                                                              // PQ multiplier has been reset.
            linearIndependentColumnsOfPQ.append(S_wavedInizializationIndex - whereQColumnChanged)
            whereQColumnChanged = S_wavedInizializationIndex
            print("\(numberOfQColumnToUse + 1)-й столбец матрицы Q будет задействован, начиная с \(S_wavedInizializationIndex + 1)-го столбца матрицы S_waved.")
        } else {
            linearIndependencyCheck += 1
            S_wavedInizializationIndex += 1
        }
    }
    print("S_waved:")
    S_waved.printBody(withPrecision: nil)
    
    print("S_waved^(-1):")
    S_waved.inversed().printBody(withPrecision: nil)
    
    print("S_waved^(-1) * P:")
    (S_waved.inversed() * P).printBody(withPrecision: nil)
    
    print("Теперь можем найти P_waved.")
    print("S_waved^(-1) * P * S_waved = P_waved:")
    let P_waved = S_waved.inversed() * P * S_waved
    P_waved.printBody(withPrecision: nil)
    
    print("Построим Q_waved.")
    print("S_waved^(-1) * Q = Q_waved:")
    let Q_waved = S_waved.inversed() * Q
    Q_waved.printBody(withPrecision: nil)
    
    print("P_waved разбивается на подсистемы c соответствующими характер. полиномами:")
    var subsystemsOfP_waved: [Matrix] = []
    var polynomialsOfSubsystemsOfP_waved: [[Double]] = []
    var upperLeftIndex = 0                                                      // Here and so forth these are used to indicate edge elements of diagonal submatrices
    var lowerRightIndex = -1
    for i in 0..<linearIndependentColumnsOfPQ.count {
        upperLeftIndex = lowerRightIndex + 1
        lowerRightIndex = upperLeftIndex + linearIndependentColumnsOfPQ[i] - 1
        subsystemsOfP_waved.append(P_waved.submatrix(from: (upperLeftIndex, upperLeftIndex), to: (lowerRightIndex, lowerRightIndex)))
        print("P\(i)\(i)_waved:")
        subsystemsOfP_waved[i].printBody(withPrecision: nil)
        polynomialsOfSubsystemsOfP_waved.append(subsystemsOfP_waved[i].characteristicPolynomial)
        print("Хар. полином P\(i)\(i)_waved:")
        printPolynomial(polynomialsOfSubsystemsOfP_waved[i], rightSideOfEquation: nil, singleVariable: true)
    }
    print("Теперь построим матрицы T для каждой из подсистем:")
    var T = Matrix(dimension: P.dimension, body: nil)                           // Consists of T matrices built from ch. polynomials of P_waved's subsystems
    upperLeftIndex = 0
    lowerRightIndex = -1
    for subsystemIndex in 0..<subsystemsOfP_waved.count {
        upperLeftIndex = lowerRightIndex + 1
        lowerRightIndex = upperLeftIndex + linearIndependentColumnsOfPQ[subsystemIndex] - 1
        for i in 0..<T.dimension.y {
            for j in 0..<T.dimension.x {
                if i == j {
                    T[i, j] = 1
                } else if i < j && i >= upperLeftIndex && i <= lowerRightIndex && j >= upperLeftIndex && j <= lowerRightIndex {
                    T[i, j] = polynomialsOfSubsystemsOfP_waved[subsystemIndex][j - i]   // [i, j] is between two edge elements of the subsystem
                }
            }
        }
        print("T\(subsystemIndex)\(subsystemIndex):")
        T.submatrix(from: (upperLeftIndex, upperLeftIndex), to: (lowerRightIndex, lowerRightIndex)).printBody(withPrecision: nil)
    }
    print("Получившаяся матрица T:")
    T.printBody(withPrecision: nil)
    
    var targetPolynomial: [Double] = [1, -targetEigenvalues[0]]                 // This polynomial will be used during a check of the answer
    for i in 1..<targetEigenvalues.count {
        targetPolynomial = polynomialsProduct(targetPolynomial, [1, -targetEigenvalues[i]])
    }
    
    print("Теперь построим матрицу Г. Для этого для каждой подсистемы вычислим векторы-строки разностей их хар. полиномов и соответствующего эталонного многочлена:")
    var G = Matrix(dimension: (P.dimension.y, subsystemsOfP_waved.count), body: nil)
    var g: Matrix                                                               // Temporary matrix that stores subsystems of G for further printing
    
    var targetSubEigenvalues: [Double]                                          // Subsystem's eigenvalues
    
    var targetSubPolynomial: [Double]                                           // Target polynomial that is built from subsystem's eigenvalues
    
    var diagonalElementInG = 0                                                  // Start point for writing information about a particular subsystem into G
    for j in 0..<G.dimension.x {
        
        g = Matrix(dimension: (1, polynomialsOfSubsystemsOfP_waved[j].count - 1), body: nil)
        
        targetSubEigenvalues = []
        for elementToAdd in diagonalElementInG..<(diagonalElementInG + polynomialsOfSubsystemsOfP_waved[j].count - 1) {
            targetSubEigenvalues.append(targetEigenvalues[elementToAdd])
        }
        targetSubPolynomial = [1, -targetSubEigenvalues[0]]
        for i in 1..<targetSubEigenvalues.count {
            targetSubPolynomial = polynomialsProduct(targetSubPolynomial, [1, -targetSubEigenvalues[i]])
        }
        
        print("Хар. полином \(j)-й подсистемы:")
        printPolynomial(polynomialsOfSubsystemsOfP_waved[j], rightSideOfEquation: nil, singleVariable: true)
        print("Эталонный многочлен \(j)-й подсистемы:")
        printPolynomial(targetSubPolynomial, rightSideOfEquation: nil, singleVariable: true)
        
        for i in diagonalElementInG..<(diagonalElementInG + polynomialsOfSubsystemsOfP_waved[j].count - 1) {
            G[i, j] = polynomialsOfSubsystemsOfP_waved[j][i - diagonalElementInG + 1] - targetSubPolynomial[i - diagonalElementInG + 1]
            g[0, i - diagonalElementInG] = G[i, j]
        }
        print("Вектор разности полиномов \(j)-й подсистемы:")
        g.printBody(withPrecision: nil)
        
        diagonalElementInG = polynomialsOfSubsystemsOfP_waved[j].count - 1
    }
    print("Расставив векторы разности по диагонали, получим матрицу Г:")
    G.printBody(withPrecision: nil)
    print("Осталось только найти матрицу М. Её формула: M = Г' * (S_waved * T)^(-1).\n")
    print("S_waved * T:")
    (S_waved * T).printBody(withPrecision: nil)
    print("(S_waved * T)^(-1):")
    (S_waved * T).inversed().printBody(withPrecision: nil)
    print("M:")
    M = G.transposed() * (S_waved * T).inversed()
    M.printBody(withPrecision: nil)
    print("Проверим найденное управление. Для этого вычислим характеристический полином матрицы P + (Q * M) и сравним его корни с корнями эталонного полинома.\n")
    print("Q * M:")
    (Q * M).printBody(withPrecision: nil)
    print("P + (Q * M):")
    (P + Q * M).printBody(withPrecision: nil)
    print("Характеристический полиной этой матрицы:")
    printPolynomial((P + Q * M).characteristicPolynomial, rightSideOfEquation: nil, singleVariable: true)
    print("Эталонный полином:")
    printPolynomial(targetPolynomial, rightSideOfEquation: nil, singleVariable: true)
    if targetPolynomial == (P + Q * M).characteristicPolynomial {
        print("Они совпадают! Все вычисления произведены верно.")
    } else {
        print("Они не совпадают! А должны. Убедитесь, что это не ошибка погрешности вычислений, и попробуйте ввести начальные данные снова.")
    }
    
}

func partiallyControllableWithScalarInput(_ P: Matrix, _ Q: Matrix, _ targetEigenvalues: [Double]) {
    print("Решается задача стабилизации системы со СКАЛЯРНЫМ входом для случая НЕПОЛНОЙ управляемости.\n")
    
    var M: Matrix
    print("Построим матрицу S:")
    var S = Matrix(dimension: P.dimension, body: nil)
    var PQ = Q
    for i in 0..<P.dimension.x {
        S.setColumn(at: i, PQ.columns[0])
        PQ = P * PQ
    }
    S.printBody(withPrecision: nil)
    let k = S.rank
    var rankS = k
    print("Проверим систему на полную управляемость: rank S = \(k) < n. Система не полностью управляема.")
    
    print("Добираем необходимое для лин. независимости количество ортов - получаем матрицу T:")
    for i in 0..<P.dimension.y {
        var e: [Double] = Array(repeating: 0, count: P.dimension.y)
        e[i] = 1
        S.setColumn(at: rankS, e)
        if rankS < S.rank {
            rankS = S.rank
            if rankS == P.dimension.y {
                break
            }
        }
    }
    let T: Matrix = S
    T.printBody(withPrecision: nil)
    
    print("T^(-1):")
    T.inversed().printBody(withPrecision: nil)
    print("P * T:")
    (P * T).printBody(withPrecision: nil)
    print("Вычислим матрицы P_waved и Q_waved.")
    print("T^(-1) * P * T = P_waved:")
    let P_waved = T.inversed() * P * T
    P_waved.printBody(withPrecision: nil)
    print("T^(-1) * Q = Q_waved:")
    let Q_waved = T.inversed() * Q
    Q_waved.printBody(withPrecision: nil)
    
    print("P_waved разбивается на подсистемы. Проверим неуправляемую подсистему на устойчивость. Вычислим характеристический полином матрицы P22_waved:")
    let P22_waved = P_waved.submatrix(from: (k, k), to: (P_waved.dimension.y - 1, P_waved.dimension.x - 1))
    P22_waved.printBody(withPrecision: nil)
    print("Характеристический полином этой матрицы:")
    printPolynomial(P22_waved.characteristicPolynomial, rightSideOfEquation: nil, singleVariable: true)
    if P22_waved.isStable() == false {
        print("Неуправляемая подсистема матрицы P_waved не стабилизируема. Следовательно, исходная система тоже не стабилизируема. Программа сейчас завершит работу..")
        return
    }
    print("Его корни находятся в левой комплексной полуплоскости. Подсистема стабилизируема.")
    
    print("Найдем характеристический полином управляемой подсистемы P11_waved:")
    let p11 = P_waved.submatrix(from: (0, 0), to: (k - 1, k - 1)).characteristicPolynomial
    printPolynomial(p11, rightSideOfEquation: nil, singleVariable: true)
    
    print("Из его коэффициентов построим матрицу K1:")
    var K1 = Matrix(dimension: (k, k), body: nil)
    for i in 0..<K1.dimension.y {
        for j in 0..<K1.dimension.x {
            if i == j {
                K1[i, j] = 1
            } else if i < j {
                K1[i, j] = p11[j - i]
            }
        }
    }
    K1.printBody(withPrecision: nil)
    
    print("Построим матрицу K, состоящую из K1 и K2(единичной):")
    var K = Matrix(dimension: P.dimension, body: nil)
    for i in 0..<K.dimension.y {
        for j in 0..<K.dimension.x {
            if i < k && j < k {
                K[i, j] = K1[i, j]
            } else if i == j {
                K[i, j] = 1
            }
        }
    }
    K.printBody(withPrecision: nil)
    
    print("Выпишем эталонный многочлен:")
    var targetPolynomial: [Double] = [1, -targetEigenvalues[0]]
    for i in 1..<targetEigenvalues.count {
        targetPolynomial = polynomialsProduct(targetPolynomial, [1, -targetEigenvalues[i]])
    }
    printPolynomial(targetPolynomial, rightSideOfEquation: nil, singleVariable: true)
    
    print("Теперь вычислим вектор строку разностей коеффициентов хар. полиномов Г:")
    var G = Matrix(dimension: (1, k), body: nil)
    for i in 0..<k {
        G[0, i] = p11[i + 1] - targetPolynomial[i + 1]
    }
    G.printBody(withPrecision: nil)
    let newRow = G.rows[0] + Array(repeating: 0, count: (P.dimension.x - k))
    G = Matrix(dimension: (1, newRow.count), body: newRow)
    print("Добьём Г нулями до строки длиной n:")
    G.printBody(withPrecision: nil)
    
    print("T * K:")
    (T * K).printBody(withPrecision: nil)
    print("(T * K)^(-1):")
    (T * K).inversed().printBody(withPrecision: nil)
    M = G * (T * K).inversed()
    print("Г * (T * K)^(-1) = M:")
    M.printBody(withPrecision: nil)
    print("Проверим найденное стабилизирующее управление.Сравним хар. полином матрицы P + (Q * M) и произведение эталонного хар. полинома на хар. полином неуправляемой подсистемы P22_waved.\n")
    print("Первый полином:")
    printPolynomial((P + Q * M).characteristicPolynomial, rightSideOfEquation: nil, singleVariable: true)
    print("Второй полином:")
    let check: [Double] = polynomialsProduct(P22_waved.characteristicPolynomial, targetPolynomial)
    printPolynomial(check, rightSideOfEquation: nil, singleVariable: true)
    if (P + Q * M).characteristicPolynomial == check {
        print("Молодец! Полиномы совпали. Вычисления произведены верно.")
    } else {
        print("Красавчик! Полиномы не совпали. Если дело не в погрешности округления, имеет смысл попробовать ввести начальные данные заново.")
    }
}

var butler = Butler()
