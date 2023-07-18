/**OPTIMAL CONTROL**

This program is supposed to be used in an online compiler such as onlineGDB.
 
BRIEF GUIDE:
 The program solves optimal stabilization problems using Zubov's method,
 providing detailed information about intermediate calculations.
 
 
INITIALIZATION:
 Matlab-like input style is used here, which means that elements in a single
 row may or may not be separated by commas:
    
    (#1) 1 2; 3 4; 5 6 [Enter] -->  [ 1, 2;
                                      3, 4;
    (#2) 1,2; 3,4; 5,6 [Enter] -->    5, 6 ]
 
 and after each row there should be a semicolon indicating its ending(except the last one,
 where the end sign is optional):
 
    (#3) 1 2 3 4; 5 6 7 8 [Enter]  -->  [ 1, 2, 3, 4;
    (#4) 1 2 3 4; 5 6 7 8; [Enter] -->    5, 6, 7, 8 ].
 
 Square brackets are also available, but why would you use them, huh?
 (#1) is recommended as the most succinct one.
 
 
 KNOWN ERRORS:
    1) Since M0 is found randomly, sometimes eigenvalues are very close to the right complex
half-plane. For example,
 
    Re(lambda) = -1*10^(-12).
 
This seems to cause problems with rounding error and stability of the system.
The issue is likely to be solved after restart of the program.

    2) If you type some characters like " ` " or " ' " while entering a matrix,
you will get an error even if the character will be deleted before you press Enter.
If this happens, just restart the program.
 
 
 OUTPUT:
 The result of the program is M_i'th, where 'i' varies from 1 to 2 inclusive.
 This behavior can be changed at the Butler structure's initialization code.
    The assistive messages in terminal guide through the whole process,
 though certain awareness obout Zubov's method is required. To achieve it you can use
 this manual:
 "Оптимальная стабилизация линейных систем:
 Учеб. пособие / Тамасян Г. Ш., Фоминых А. В. - СПб.:
 Изд-во ВВМ, 2022. - 66 с. Библ. 14 назв."
 ISBN 978-5-9651-1405-4
 
*/

                                                                                            // Efficient determinant; somewhat efficient inversion... A convenient matrix structure for daily usage
struct Matrix {
    var body: [Double] = []
    var dimension: (y: Int, x: Int)
    
                                                                                            // Returns the whole matrix in the format of a 2-dimensional array
    var rows: [[Double]] {
        var answer: [[Double]] = Array(repeating: Array(repeating: 0, count: dimension.x), count: dimension.y)
        for i in 0..<dimension.y {
            for j in 0..<dimension.x {
                answer[i][j] = self[i, j]
            }
        }
        return answer
    }
    
                                                                                            // The same as rows, but columns
    var columns: [[Double]] {
        var answer: [[Double]] = Array(repeating: Array(repeating: 0, count: dimension.y), count: dimension.x)
        for j in 0..<dimension.x {
            for i in 0..<dimension.y {
                answer[j][i] = self[i, j]
            }
        }
        return answer
    }
    
                                                                                            // Uses Danilevsky method to obtain coefficients of the characteristic polynomial
    var characteristicPolynomial: [Double] {
        var B = Matrix(dimension: self.dimension, body: nil)
        var D = Matrix(dimension: self.dimension, body: self.body)                          // A matrix to be transformed into Frobenius normal form
        var answer: [Double] = []
        let n = self.dimension.y - 1
        for step in 1...(self.dimension.y - 1) {
                                                                                            //k = (n - step + 1) - a transition function between what I'd read and what I wrote
            
            if D[n - step + 1, n - step] == 0 {                                             // The element is going to be used as a divider. This is a try to exchange it
                for m in 0..<(n - step) {                                                   // for another one at the left(with smaller second index) using determinant rules
                    if D[n - step + 1, m] != 0 {
                        D.swapColumns(n - step, m)
                        D.swapRows(n - step, m)
                        break
                    }
                }
                                                                                            // If there is no available elements for an exchange, an original
                if D[n - step + 1, n - step] == 0 {                                         // polynomial can be calculated as a product of two:
                    var det1: [Double]                                                      // Upper left...
                    if D.submatrix(from: (0, 0), to: (n - step, n - step)).body.count == 1 {
                        det1 = [1, -D[0, 0]]
                    } else {
                        det1 = D.submatrix(from: (0, 0), to: (n - step, n - step)).characteristicPolynomial
                    }
                    var det2: [Double] = [1]                                                // ...and lower right
                    for k in (n - step + 1)..<D.dimension.x {
                        det2.append(-D[(n - step + 1), k])
                    }
                    var det1det2: [Double] = Array(repeating: 0, count: det1.count + det2.count - 1) // Product of two polynomials
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
            for i in 0..<B.dimension.y {                                                    // A normal process with non-zero divider on the way of the Danilevsky method
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
    
    init(dimension: (y: Int, x: Int), body: [Double]?) {
        if let body = body {
            self.body = body
        } else {
            self.body = Array(repeating: 0, count: dimension.x * dimension.y)
        }
        self.dimension = dimension
    }
    
                                                                                            // Removes [], if there are any, split the string by ';' into rows
                                                                                            // and split rows into elements by ' ', ',' and ', '
    init(sBody: String) {
        var bodyWithoutSquareBraces = sBody
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
    
    subscript(_ y: Int, _ x: Int) -> Double {
        get {
            return body[y * dimension.x + x]
        }
        set {
            body[y * dimension.x + x] = newValue
        }
    }
    
                                                                                            // Matrix + Matrix
    static func + (mLeft: Matrix, mRight: Matrix) -> Matrix {
        var answer = Matrix(dimension: mLeft.dimension, body: nil)
        for i in 0..<answer.body.count {
            answer.body[i] = mLeft.body[i] + mRight.body[i]
        }
        return answer
    }
    
                                                                                            // Matrix - Matrix
    static func - (mLeft: Matrix, mRight: Matrix) -> Matrix {
        var answer = Matrix(dimension: mLeft.dimension, body: nil)
        for i in 0..<answer.body.count {
            answer.body[i] = mLeft.body[i] - mRight.body[i]
        }
        return answer
    }
    
                                                                                            // Double * Matrix
    static func * (dLeft: Double, mRight: Matrix) -> Matrix {
        var answer = Matrix(dimension: mRight.dimension, body: nil)
        for i in 0..<mRight.body.count {
            answer.body[i] = dLeft * mRight.body[i]
        }
        return answer
    }
    
                                                                                            // Matrix * Matrix
    static func * (mLeft: Matrix, mRight: Matrix) -> Matrix {
        if mLeft.body.count == 1 {                                                          // If one of the operands is just a number(a single-element matrix) ->
            return mLeft[0, 0] * mRight                                                     // delegate the call to an appropriate function
        }
        if mRight.body.count == 1 {
            return mRight[0, 0] * mLeft
        }
        
        var answer = Matrix(dimension: (mLeft.dimension.y, mRight.dimension.x), body: nil)
        for i in 0..<answer.dimension.y {                                                   // Real matrix product
            for j in 0..<answer.dimension.x {
                for k in 0..<mLeft.dimension.x {
                    answer[i, j] += mLeft[i, k] * mRight[k, j]
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
    
                                                                                            // Returns submatrix which borders in original matrix
                                                                                            // are defined by two elements: upper left and lower right
    func submatrix(from upperLeft: (y: Int, x: Int), to lowerRight: (y: Int, x: Int)) -> Matrix {
        var answer = Matrix(dimension: (lowerRight.y - upperLeft.y + 1, lowerRight.x - upperLeft.x + 1), body: nil)
        for i in upperLeft.y...lowerRight.y {
            for j in upperLeft.x...lowerRight.x {
                answer[i - upperLeft.y, j - upperLeft.x] = self[i, j]
            }
        }
        return answer
    }
    
                                                                                            // For 1x1, 2x2 and 3x3 is calculated by explicit formulas.
                                                                                            // 4x4+ is calculated with modification of the Gauss method
    func det() -> Double {
        if dimension.x == 1 && dimension.y == 1 {
            return self[0, 0]
        } else if dimension.x == 2 && dimension.y == 2 {
            return self[0, 0] * self[1, 1] - self[0, 1] * self[1, 0]
        } else if dimension.x == 3 && dimension.y == 3 {
            return self[0, 0] * self[1, 1] * self[2, 2] + self[1, 0] * self[2, 1] * self[0, 2] + self[2, 0] * self[0, 1] * self[1, 2] - self[2, 0] * self[1, 1] * self[0, 2] - self[1, 0] * self[0, 1] * self[2, 2] - self[0, 0] * self[2, 1] * self[1, 2]
        } else {
            var temporaryMatrix = self
            var answer: Double = 1
            for d in 0..<temporaryMatrix.dimension.x {                                      // Excludes zero elements from the diagonal
                if temporaryMatrix[d, d] == 0 {
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
                for i in (d+1)..<temporaryMatrix.dimension.y {                              // Forward elimination from the Gauss method
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
    
                                                                                            // 1x1, 2x2 and 3x3 are computed by explicit formulas.
                                                                                            // 4x4+ inversion is based on the Cramer's rule
    func inversed() -> Matrix {
        var answer = Matrix(dimension: dimension, body: nil)
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
                    var submatrixBody: [Double] = []                                        // a minor of a smaller dimension
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
    
                                                                                            // Uses Routh-Hurwitz stability criterion and characteristicPolynomial propety of the matrix structure
    func isStable() -> Bool {
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
    
                                                                                            // Note: rounding is not mathematical: it just drops all signs after the specified one
    func printBody(withPrecision: Int?) {
        for i in 0..<dimension.y {
            for j in 0..<dimension.x {
                var element = String(self[i, j])
                if var precision = withPrecision {
                    if precision == 0 {                                                     // This erases a dot between integer and fractional parts of a number
                        precision = -1
                    }
                    
                    if element.contains(where: {$0 == "e" || $0 == "E"}) {                  // Processing numbers in scientific format
                        var temporaryString = element
                        
                        temporaryString.remove(at: temporaryString.firstIndex(of: ".")!)    // Extract non-zero digits from the number
                        temporaryString = String(temporaryString[..<temporaryString.firstIndex(of: "e")!])
                        
                                                                                            // Combine the string from the previous step with the proper amount of
                                                                                            // zeros, acquired from an exponent of the original number
                        element = "0." + String(repeating: "0", count: Int(element[element.index(element.firstIndex(of: "e")!, offsetBy: 2)...])! - 1) + temporaryString
                    }
                    
                    let indexOfTheDot = element.firstIndex(of: ".")                         // Precision adjustments
                    
                    if precision - element[element.firstIndex(of: ".")!...].count >= 0 {    // If string representation of the elements is shorter than
                        element += String(repeating: "0", count: precision - element[element.firstIndex(of: ".")!...].count + 1)    //required precision, it is filled with zeros
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

                                                                                            // An auxiliary structure that takes care about human-readable initialization process
                                                                                            // and some other visual amenities
struct Butler {
    var P: Matrix                                                                           // Matrices that are known from the start
    var Q: Matrix
    var A: Matrix
    var B: Matrix
    var C: Matrix
    var M: Matrix
    
    let n: Int = 2                                                                          // Number of iterations to run
    
    init() {
        var input: String
        print("P: ", terminator: "")
        input = readLine()!
        P = Matrix(sBody: input)
        
        print("Q: ", terminator: "")
        input = readLine()!
        Q = Matrix(sBody: input)
        
        print("A: ", terminator: "")
        input = readLine()!
        A = Matrix(sBody: input)
        
        print("B: ", terminator: "")
        input = readLine()!
        B = Matrix(sBody: input)
        if B.body == [0] {
            B = Matrix(dimension: (Q.dimension), body: nil)
        }
        
        print("C: ", terminator: "")
        input = readLine()!
        C = Matrix(sBody: input)
        
        print("M(если нач. управление неизвестно, нажмите Enter):")
        input = readLine()!
        if input != "" {
            M = Matrix(sBody: input)
        } else {
            M = someFeasibleM(P, Q)
            print("Рандомное допустимое M0:")
            M.printBody(withPrecision: 3)
            print("Матрица (P + QM0) =")
            (P + Q * M).printBody(withPrecision: 3)
            print("Её характеристический полином:")
            printPolynomial((P + Q * M).characteristicPolynomial, rightSideOfEquation: nil, singleVariable: true)
        }
        
        for i in 1...n {
            print("\(i)-АЯ ИТЕРАЦИЯ:")
            M = ZubovMethod(P, Q, A, B, C, M)
            print("C^(-1) * (Theta * Q - B)^T = M\(i) =")
            M.printBody(withPrecision: 3)
        }
        print("Проверим допустимость найденного управления:")
        print("Q * M2 =")
        (Q * M).printBody(withPrecision: 3)
        print("P + QM2 =")
        (P + Q * M).printBody(withPrecision: 3)
        print("Характеристический полином (P + QM2):")
        printPolynomial((P + Q * M).characteristicPolynomial, rightSideOfEquation: nil, singleVariable: true)
        if (P + Q * M).isStable() {
            print("Поздравляю! Система стабилизирована. Все собственные числа лежат в левой косплексной полуплоскости.")
        } else {
            print("Поздравляю! Система не стабилизирована. Это не норм. Попробуйте перезапустить программу, введя начальные данные с большей внимательностью.")
        }
    }
}

func ZubovMethod(_ P: Matrix, _ Q: Matrix, _ A: Matrix, _ B: Matrix, _ C: Matrix, _ M: Matrix) -> Matrix {
    var Theta = Matrix(dimension: (P.dimension.y, P.dimension.x), body: nil)
    print("Решим матричное уравнение Ляпунова.")
    let PQM: Matrix = (P + Q * M)
    print("QM =")
    (Q * M).printBody(withPrecision: 3)
    print("P + QM =")
    PQM.printBody(withPrecision: 3)
    let ABMBMMCM: Matrix = (A + B * M + (B * M).transposed() + M.transposed() * (C * M))
    print("BM =")
    (B * M).printBody(withPrecision: 3)
    print("(M^T)CM =")
    (M.transposed() * C * M).printBody(withPrecision: 3)
    print("A + BM + (BM)^T + (M^T)CM =")
    ABMBMMCM.printBody(withPrecision: 3)
    
                                                                                            //The problem requires to solve a system of linear equations. The solution for 2x2 and
                                                                                            // 3x3 cases can be found explicitly. 4x4 is not included
    var thetaSLE_left = Matrix(sBody: "0")
    var thetaSLE_right = Matrix(sBody: "0")                                                 // SLE stands for System of Linear Equations
    if P.dimension.y == 2 && P.dimension.x == 2 {
        thetaSLE_left = Matrix(dimension: (3, 3), body: [2 * PQM[0, 0], 2 * PQM[1, 0], 0, PQM[0, 1], PQM[1, 1] + PQM[0, 0], PQM[1, 0], 0, 2 * PQM[0, 1], 2 * PQM[1, 1]])
        thetaSLE_right = Matrix(dimension: (3, 1), body: [ABMBMMCM[0, 0], ABMBMMCM[0, 1], ABMBMMCM[1, 1]])
        let thetaElements = thetaSLE_left.inversed() * thetaSLE_right
        Theta[0, 0] = thetaElements[0, 0]
        Theta[0, 1] = thetaElements[1, 0]
        Theta[1, 0] = thetaElements[1, 0]
        Theta[1, 1] = thetaElements[2, 0]
    } else if P.dimension.y == 3 && P.dimension.x == 3 {
                                                                                            // A pretty sick line
        thetaSLE_left = Matrix(dimension: (6, 6), body: [2 * PQM[0, 0], 2 * PQM[1, 0], 2 * PQM[2, 0], 0, 0, 0, PQM[0, 1], PQM[1, 1] + PQM[0, 0], PQM[2, 1], PQM[1, 0], PQM[2, 0], 0, PQM[0, 2], PQM[1, 2], PQM[2, 2] + PQM[0, 0], 0, PQM[1, 0], PQM[2, 0], 0, 2 * PQM[0, 1], 0, 2 * PQM[1, 1], 2 * PQM[2, 1], 0, 0, PQM[0, 2], PQM[0, 1], PQM[1, 2], PQM[2, 2] + PQM[1, 1], PQM[2, 1], 0, 0, 2 * PQM[0, 2], 0, 2 * PQM[1, 2], 2 * PQM[2, 2]])
        thetaSLE_right = Matrix(dimension: (6, 1), body: [ABMBMMCM[0, 0], ABMBMMCM[0, 1], ABMBMMCM[0, 2], ABMBMMCM[1, 1], ABMBMMCM[1, 2], ABMBMMCM[2, 2]])
        let thetaElements = thetaSLE_left.inversed() * thetaSLE_right
        Theta[0, 0] = thetaElements[0, 0]
        Theta[0, 1] = thetaElements[1, 0]
        Theta[0, 2] = thetaElements[2, 0]
        Theta[1, 0] = thetaElements[1, 0]
        Theta[1, 1] = thetaElements[3, 0]
        Theta[1, 2] = thetaElements[4, 0]
        Theta[2, 0] = thetaElements[2, 0]
        Theta[2, 1] = thetaElements[4, 0]
        Theta[2, 2] = thetaElements[5, 0]
        
    } else {                                                                                // The 4x4+ matrix equation solver isn't needed for this task: you will only get 2x2 and 3x3 matrices
        print(" Алгоритм решения задачи для матрицы 4-го порядка не реализовано в этой версии Антизаранника. Программа сейчас завершит работу... ")
        return Matrix(sBody: "")                                                            // This will cause an error
    }
    print("Получили систему уравнений относительно элементов Theta:")
    for i in 0..<thetaSLE_left.dimension.y {
        printPolynomial(thetaSLE_left.rows[i], rightSideOfEquation: thetaSLE_right.rows[i], singleVariable: false)
    }
    print("Решив систему, получаем искомую матрицу Theta:")
    Theta.printBody(withPrecision: 3)
    print("Theta * Q =")
    (Theta * Q).printBody(withPrecision: 3)
    print("Theta * Q - B =")
    (Theta * Q - B).printBody(withPrecision: 3)
    let newM = C.inversed() * (Theta * Q - B).transposed()
    return newM
}

                                                                                            // Is used for writing accompanying information about calculations into terminal window
func printPolynomial(_ polynomial: [Double], rightSideOfEquation: [Double]?, singleVariable: Bool) {
    
                                                                                            // From there till the end just some visual stuff' manipulations
    
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

                                                                                            // Finds random feasible starting integer-valued M0 for Zubov's method.
func someFeasibleM(_ P: Matrix, _ Q: Matrix) -> Matrix {
    var M = Matrix(dimension: (Q.dimension.x, Q.dimension.y), body: nil)
    
    for attempt in 0...(10000 * M.dimension.x * M.dimension.y) {
        
        for i in 0..<M.dimension.y {
            for j in 0..<M.dimension.x {

                M[i, j] = Double(Int.random(in: -(attempt/2000 + 1)...(attempt/2000) ))     // Limits for random numbers to use in the matrix are increasing with
            }                                                                               // the number of unsuccessful attempts
        }
        if (P + Q * M).isStable() {
            return M
        }
    }
    print("ДОПУСТИМОГО НАЧАЛЬНОГО ПРИБЛИЖЕНИЯ НЕ СУЩЕСТВУЕТ - МАТРИЦА 'P' ИСХОДНОЙ СИСТЕМЫ НЕ МОЖЕТ БЫТЬ СТАБИЛИЗИРОВАНА. ПРОГРАММА СЕЙЧАС ЗАВЕРШИТ РАБОТУ.")   //Well, that's just a guess, but ok
    return Matrix(dimension: (0, 0), body: nil)
}

var butler = Butler()



