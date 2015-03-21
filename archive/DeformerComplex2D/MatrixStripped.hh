////
// MatrixStripped.hh
////
// Author: Julian Panetta
// Description:
//		Describes a matrix class with addition, subtraction, multiplication,
//		and inverstion operations.
// Revision History:
//		10/23/07 - Initial Revision
//		10/27/07 - Added overloaded operators
//                 removed methods the operators replaced
//      10/30/07 - Changed bounds checking on element accessors/mutators
//                 to restrict access to undefined elements, not just
//                 unallocated memory addresses. Minor style changes
//      10/08/08 - Generalized into a template for use with other data types
//                 Added function operator overload to implement indexing.
//                 Added scalar multiplication, and inversion.
//      10/14/08 - Added transposition, normalization, and homogenization
//      10/27/08 - Made normalization only normalize the first n-1 elemnts
//                 of an n-vector (our vectors are [x, y, z, w] and we only
//                 want to normalize [x, y, z])
//      01/31/09 - Modified to work with the Vector2 complex number class
////

#ifndef MATRIX_LIBRARY_INCLUDED
#define MATRIX_LIBRARY_INCLUDED

#include <cassert>
#include <math.h>
#include <stdio.h>

#define tolerance 1.0e-4

template<typename T> class MatrixStripped
{
    private:
        int rows, cols;             // The number of rows and columns
        T *data;                    // The matrix data
        
        // Helper functions
        void copy(const MatrixStripped<T> &aMatrix);
        void cleanup();

        bool zero(T value) const;         // Check if a value is zero

    public:
        // Constructors
        MatrixStripped()    {                       // Default constructor (0x0 matrix)
            rows = cols = 0;
            data = new T[0];
        }

        MatrixStripped(int r, int c)    {           // New matrix constructor
            assert( (r >= 0) && (c >= 0));

            rows = r;
            cols = c;
            data = new T[r*c];
            
            /* NOTE: I AM NO LONGER INITIALIZING NEWLY CONSTRUCTED MATRICES TO ZERO
             * BECAUSE THE I CANNOT ASSIGN 0 TO A VECTOR2 */
            //for (int i = 0; i < r * c; i++)
            //  data[i] = 0;
        }

        MatrixStripped(const MatrixStripped<T> &aMatrix)   {  // Copy constructor
            copy(aMatrix);
        }
        
        ~MatrixStripped()   {   // Destructor
            cleanup();
        }
        
        // Mutator Methods
        void setelem(int r, int c, T value);
        void swapRows(int r1, int r2);
        void subtractMultiple(int r1, int r2, T scalar);
        void multiplyRow(int r, T scalar);

        void normalize();
        void homogenize();

        // Accessor Methods
        int getrows() const;
        int getcols() const;
        T getelem(int r, int c) const;
        const MatrixStripped<T> transpose() const;

        // Print Method
        void print() const;
        
        // Overloaded Operators
        T &operator()(int r, int c);                        // Function operator
        const T &operator()(int r, int c) const;
        MatrixStripped<T> & operator=(const MatrixStripped<T> &aMatrix);    // Assignment operator
        MatrixStripped<T> & operator*=(const MatrixStripped<T> &aMatrix);
        MatrixStripped<T> & operator*=(const T scalar);
        MatrixStripped<T> & operator+=(const MatrixStripped<T> &aMatrix);
        MatrixStripped<T> & operator-=(const MatrixStripped<T> &aMatrix);

        const MatrixStripped<T> operator*(const MatrixStripped<T> &aMatrix) const;
        const MatrixStripped<T> operator*(const T scalar) const;
        const MatrixStripped<T> operator+(const MatrixStripped<T> &aMatrix) const;
        const MatrixStripped<T> operator-(const MatrixStripped<T> &aMatrix) const;

        // const MatrixStripped<T> operator*(const T scalar) const;

        bool operator==(const MatrixStripped<T> &aMatrix) const;
        bool operator!=(const MatrixStripped<T> &aMatrix) const;
};

////
// Print matrix to stdout
////
template<typename T>
void MatrixStripped<T>::print() const
{
    for (int i = 0; i < getrows(); i++)
        for (int j = 0; j < getcols(); j++)
            printf("%s%f%s", (j == 0) ? "[ " : "", (double)getelem(i, j),
                    (j != getcols() - 1) ? ", " : " ]\n");
}

////
// Mutator Methods
////

// Set the integer value of entry Matrix[r,c]
template<typename T>
void MatrixStripped<T>::setelem(int r, int c, T value)
{
	// Prevent write to undefined elements.	
	assert((r >= 0 ) && (c >= 0) && (r < rows) && (c < cols));
	
	data[r * cols + c] = value;
}

////
// Accessor Methods
////

template<typename T>
int MatrixStripped<T>::getrows() const
{
	return rows;
}

template<typename T>
int MatrixStripped<T>::getcols() const
{
	return cols;
}

// Retrieve the integer at row r and column c.
template<typename T>
T MatrixStripped<T>::getelem(int r, int c) const
{
	// Prevent access to undefined elements.
	assert((r >= 0 ) && (c >= 0) && (r < rows) && (c < cols));
	
	return data[r * cols + c];
}

////
// Helper Functions
////
template<typename T>
void MatrixStripped<T>::copy(const MatrixStripped<T> &aMatrix)
{
	// Allocate array of the same dimensions as aMatrix.
	rows = aMatrix.rows;
	cols = aMatrix.cols;
	data = new T[rows * cols];

	// Copy the elmeents from aMatrix into the new array.
	for (int i = 0; i < rows * cols; i++)
		data[i] = aMatrix.data[i];
}

template<typename T>
void MatrixStripped<T>::cleanup()
{
	// The data array was dynamically allocated and it must be deallocated
	delete[] data;
}

// Decides whether an element is zero using a given tolerance
template<typename T>
bool MatrixStripped<T>::zero(T value) const
{
    return (fabs(value) < tolerance);
}

////
// Overloaded Operators
////

// Function call operator overload (constant reference)
template<typename T>
const T &MatrixStripped<T>::operator()(int r, int c) const
{
	// Prevent access to undefined elements.
	assert((r >= 0 ) && (c >= 0) && (r < rows) && (c < cols));

    return data[r * cols + c];
}

// Function call operator overload (reference)
template<typename T>
T &MatrixStripped<T>::operator()(int r, int c)
{
	// Prevent access to undefined elements.
	assert((r >= 0 ) && (c >= 0) && (r < rows) && (c < cols));

    return data[r * cols + c];
}

// Assignment operator
template<typename T>
MatrixStripped<T> & MatrixStripped<T>::operator=(const MatrixStripped<T> &aMatrix)
{
	// Prevent self-assignment
	if (this != &aMatrix)	{
		cleanup();
		copy(aMatrix);
	}
	return *this;
}

// Matrix multiplication operator
//	Multiplication code is held within this simple operator so that an
//	extra copy need not be made:
//	-Currently one copy of *this is made to accumulate multiplication results
//	-If =* held the code, then for *, a copy of *this to a temporary
//	 matrix would be required, and then a second copy to an accumulator
//	 matrix would be needed.
template<typename T>
const MatrixStripped<T> MatrixStripped<T>::operator*(const MatrixStripped<T> &aMatrix) const
{
	// Matrix multiplication A * B is only defined for cases where the number
	// of rows in A is equal to the number of columns in B.
	assert(cols == aMatrix.rows);
	// Create a temporary matrix into which we will perform multiplication
	// The dimensions are this.rows by aMatrix.cols
	int resultRows = rows, resultCols = aMatrix.cols;
	MatrixStripped resultMatrix(resultRows, resultCols);
	
	// Loop over the rows of the temporary matrix
	for (int i = 0; i < resultRows; i++)	{
		// Loop over the columns of the temporary matrix
		for (int j = 0; j < resultCols; j++)	{
			T dotProduct = Vector2(0, 0);
			// Compute value of the i,jth entry of the temporary matrix by 
			// computing the dot product of the ith row vector of this matrix
			// with the jth column vector of aMatrix
			for (int element = 0; element < cols; element++)
				dotProduct += getelem(i, element) * aMatrix.getelem(element, j);
			resultMatrix.setelem(i, j, dotProduct);
		}
	}
	
	return resultMatrix;
}

// Matrix multiply compound assignment operator
template<typename T>
MatrixStripped<T> & MatrixStripped<T>::operator*=(const MatrixStripped<T> &aMatrix)
{
	// Reuse the * operator defined above;
	return *this = *this * aMatrix;
}

// Scalar multiplication operator
template<typename T>
const MatrixStripped<T> MatrixStripped<T>::operator*(const T scalar) const
{
    MatrixStripped resultMatrix(rows, cols);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            resultMatrix.setelem(i, j, getelem(i, j) * scalar);

    return resultMatrix;
}

// Scalar multiply compound assignment operator
template<typename T>
MatrixStripped<T> & MatrixStripped<T>::operator*=(const T scalar)
{
	// Reuse the * operator defined above;
	return *this = *this * scalar;
}

// Addition compound assignment operator
template<typename T>
MatrixStripped<T> & MatrixStripped<T>::operator+=(const MatrixStripped<T> &aMatrix)
{
	// Matrix addition A + B is only defined for A, B of equal dimensions
	assert((rows == aMatrix.rows) && (cols == aMatrix.cols));

	// Add all corresponding entries of aMatrix to this matrix
	for (int i = 0; i < rows * cols; i++)
		data[i] += aMatrix.data[i];
		
	return *this;	//return this matrix to allow operator chaining
}

// Subtraction compound assignment operator
template<typename T>
MatrixStripped<T> & MatrixStripped<T>::operator-=(const MatrixStripped<T> &aMatrix)
{
	// Matrix subtraction A - B is only defined for A, B of equal dimensions
	assert((rows == aMatrix.rows) && (cols == aMatrix.cols));

	// Subtract all corresponding entries of aMatrix from this matrix
	for (int i = 0; i < rows * cols; i++)
		data[i] -= aMatrix.data[i];
	
	return *this;	//return this matrix to allow operator chaining
}

// Addition operator
template<typename T>
const MatrixStripped<T> MatrixStripped<T>::operator+(const MatrixStripped<T> &aMatrix) const
{
	// Add aMatrix to a copy of this matrix and return the result
	return MatrixStripped(*this) += aMatrix;	
}

// Subtraction operator
template<typename T>
const MatrixStripped<T> MatrixStripped<T>::operator-(const MatrixStripped<T> &aMatrix) const
{
	// Subtract aMatrix from a copy of this matrix and return the result
	return MatrixStripped(*this) -= aMatrix;
}

// Equality operator
template<typename T>
bool MatrixStripped<T>::operator==(const MatrixStripped<T> &aMatrix) const
{
	// If the dimensions do not equal, the matrices cannot equal
	if ((rows != aMatrix.rows) || (cols != aMatrix.cols))
		return false;

	// All elements in the arrays must equal for the matrices to equal.
	for (int i = 0; i < rows * cols; i++)	{
		if (data[i] != aMatrix.data[i])
			return false;
	}
	return true;
}

// Inequality operator
template<typename T>
bool MatrixStripped<T>::operator!=(const MatrixStripped<T> &aMatrix) const
{
	//reuse the equality operator
	return !(*this == aMatrix);
}

// Returns the transpose of this matrix
template<typename T>
const MatrixStripped<T> MatrixStripped<T>::transpose() const
{
    MatrixStripped transposed(cols, rows);

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            transposed(j, i) = getelem(i, j);

    return transposed;
}

#endif
