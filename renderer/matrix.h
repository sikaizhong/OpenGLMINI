#pragma once
#include<vector.h>
namespace CGL {
    using index_t = unsigned;
    template <class Scalar >
    inline Scalar det2x2(Scalar a11, Scalar a12, Scalar a21, Scalar a22) {
        return a11 * a22 - a12 * a21;
    }

    template <class Scalar >
    inline Scalar det3x3(
        Scalar a11, Scalar a12, Scalar a13,
        Scalar a21, Scalar a22, Scalar a23,
        Scalar a31, Scalar a32, Scalar a33) {
        return
            a11 * det2x2(a22, a23, a32, a33)
            - a21 * det2x2(a12, a13, a32, a33)
            + a31 * det2x2(a12, a13, a22, a23);
    }
    template <class Scalar >

    inline Scalar det4x4(
        Scalar a11, Scalar a12, Scalar a13, Scalar a14,
        Scalar a21, Scalar a22, Scalar a23, Scalar a24,
        Scalar a31, Scalar a32, Scalar a33, Scalar a34,
        Scalar a41, Scalar a42, Scalar a43, Scalar a44
    ) {
        Scalar m12 = a21 * a12 - a11 * a22;
        Scalar m13 = a31 * a12 - a11 * a32;
        Scalar m14 = a41 * a12 - a11 * a42;
        Scalar m23 = a31 * a22 - a21 * a32;
        Scalar m24 = a41 * a22 - a21 * a42;
        Scalar m34 = a41 * a32 - a31 * a42;

        Scalar m123 = m23 * a13 - m13 * a23 + m12 * a33;
        Scalar m124 = m24 * a13 - m14 * a23 + m12 * a43;
        Scalar m134 = m34 * a13 - m14 * a33 + m13 * a43;
        Scalar m234 = m34 * a23 - m24 * a33 + m23 * a43;

        return (m234 * a14 - m134 * a24 + m124 * a34 - m123 * a44);
    }

    template <unsigned DIM, class FT>
    class Matrix {
    public:
        typedef Matrix<DIM, FT> matrix_type;
        typedef FT value_type;
        static const index_t dim = DIM;
        inline Matrix() {
            load_identity();
        }
        explicit Matrix(const FT* vals) {
            for (index_t i = 0; i < DIM; i++) {
                for (index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] = *vals;
                    ++vals;
                }
            }
        }

        /**
         * \brief Gets the matrix dimension
         * \return the value of \p DIM
         */
        inline index_t dimension() const {
            return DIM;
        }

        /**
         * \brief Clears the matrix
         * \details This resets all values to 0 (zero)
         */
        inline void load_zero() {
            for (index_t i = 0; i < DIM; i++) {
                for (index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] = FT(0);
                }
            }
        }

        /**
         * \brief Sets the matrix to identity
         * \details This sets all coefficients of this matrix to be equal to
         * \p DIM x \p DIM identity matrix.
         */
        inline void load_identity() {
            for (index_t i = 0; i < DIM; i++) {
                for (index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] = (i == j) ? FT(1) : FT(0);
                }
            }
        }

        /**
         * \brief Tests whether a matrix is the identity matrix.
         * \retval true if the matrix is the identity matrix
         * \retval false otherwise
         */
        inline bool is_identity() const {
            for (index_t i = 0; i < DIM; i++) {
                for (index_t j = 0; j < DIM; j++) {
                    FT rhs = ((i == j) ? FT(1) : FT(0));
                    if (coeff_[i][j] != rhs) {
                        return false;
                    }
                }
            }
            return true;
        }

        /**
         * \brief Gets a modifiable element
         * \details Gets element at row \p i and column \p j in the matrix. If
         * indices are out of range, the function calls abort().
         * \param[in] i row index of the element
         * \param[in] j column index of the element
         * \return a reference to the element at coordinates (\p i, \p j).
         */
        inline FT& operator() (index_t i, index_t j) {
            return coeff_[i][j];
        }

        /**
         * \brief Gets a non-modifiable element
         * \details Gets element at row \p i and column \p j in the matrix. If
         * indices are out of range, the function calls abort().
         * \param[in] i row index of the element
         * \param[in] j column index of the element
         * \return a const reference to the element at coordinates (\p i, \p j).
         */
        inline const FT& operator() (index_t i, index_t j) const {
            return coeff_[i][j];
        }
        inline FT* operator[](index_t index) {
            return coeff_[index];
        }

        // Overloaded const version for read-only access
        inline const FT* operator[](index_t index) const {
            return coeff_[index];
        }
        /**
         * \brief Adds a matrix in place
         * \details This adds matrix \p m to this matrix in place.
         * \param[in] m a matrix of the same dimension
         * \return a reference to this matrix
         */
        inline matrix_type& operator+= (const matrix_type& m) {
            for (index_t i = 0; i < DIM; i++) {
                for (index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] += m.coeff_[i][j];
                }
            }
            return *this;
        }

        /**
         * \brief Subtracts a matrix in place
         * \details This subtracts matrix \p m from this matrix in place.
         * \param[in] m a matrix of the same dimension
         * \return a reference to this matrix
         */
        inline matrix_type& operator-= (const matrix_type& m) {
            for (index_t i = 0; i < DIM; i++) {
                for (index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] -= m.coeff_[i][j];
                }
            }
            return *this;
        }

        /**
         * \brief Multiplies by a scalar in place
         * \details This multiplies all the coefficients of this matrix by the
         * value \p val.
         * \param[in] val a scalar value of the same type than matrix elements
         * \return a reference to this matrix
         */
        inline matrix_type& operator*= (FT val) {
            for (index_t i = 0; i < DIM; i++) {
                for (index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] *= val;
                }
            }
            return *this;
        }

        /**
         * \brief Divides by a scalar in place
         * \details This divides all the coefficients of this matrix by the
         * value \p val.
         * \param[in] val a scalar value of the same type than matrix elements
         * \return a reference to this matrix
         */
        inline matrix_type& operator/= (FT val) {
            for (index_t i = 0; i < DIM; i++) {
                for (index_t j = 0; j < DIM; j++) {
                    coeff_[i][j] /= val;
                }
            }
            return *this;
        }

        /**
         * \brief Adds 2 matrices
         * \details Builds a matrix by adding matrix \p m to this matrix.
         * \param[in] m another matrix
         * \return the matrix (\p this + \p m)
         */
        inline matrix_type operator+ (const matrix_type& m) const {
            matrix_type result = *this;
            result += m;
            return result;
        }

        /**
         * \brief Subtracts 2 matrices
         * \details Builds a matrix by subtracting matrix \p m to this matrix.
         * \param[in] m another matrix
         * \return the matrix (\p this + \p m)
         */
        inline matrix_type operator- (const matrix_type& m) const {
            matrix_type result = *this;
            result -= m;
            return result;
        }

        /**
         * \brief Multiplies a matrix by a scalar
         * \details Builds a matrix by multiplying all the coefficients of
         * this matrix by scalar value \p val.
         * \param[in] val a scalar value of the same type than matrix elements
         * \return the resulting matrix
         */
        inline matrix_type operator* (FT val) const {
            matrix_type result = *this;
            result *= val;
            return result;
        }

        /**
         * \brief Divides a matrix by a scalar
         * \details Builds a matrix by dividing all the coefficients of
         * this matrix by scalar value \p val.
         * \param[in] val a scalar value of the same type than matrix elements
         * \return the resulting matrix
         */
        inline matrix_type operator/ (FT val) const {
            matrix_type result = *this;
            result /= val;
            return result;
        }

        /**
         * \brief Multiplies 2 matrices
         * \details Builds a matrix by multiplying this matrix by matrix \p
         * m.
         * \param[in] m another matrix
         * \return the matrix (\p this * \p m)
         */
        matrix_type operator* (const matrix_type& m) const {
            matrix_type result;
            for (index_t i = 0; i < DIM; i++) {
                for (index_t j = 0; j < DIM; j++) {
                    result.coeff_[i][j] = FT(0);
                    for (index_t k = 0; k < DIM; k++) {
                        result.coeff_[i][j] += coeff_[i][k] * m.coeff_[k][j];
                    }
                }
            }
            return result;
        }

        /**
         * \brief Computes the inverse matrix
         * \details Computes matrix \p M such that (\p this * \p M) = identity
         * \return the inverse matrix
         */
        matrix_type inverse() const {
            matrix_type result;
            bool invertible = compute_inverse(result);
            geo_assert(invertible);
            return result;
        }


        /**
         * \brief Computes the inverse matrix
         * \details Computes matrix \p M such that (\p this * \p M) = identity
         * \param[out] result the inverse matrix
         * \return true if the matrix is inversible
         * \retval false otherwise
         */
        bool compute_inverse(matrix_type& result) const {
            FT val = FT(0.0), val2 = FT(0.0);
            matrix_type tmp = (*this);

            result.load_identity();

            for (index_t i = 0; i != DIM; i++) {
                val = tmp(i, i);                     /* find pivot */
                index_t ind = i;
                for (index_t j = i + 1; j != DIM; j++) {
                    if (fabs(tmp(j, i)) > fabs(val)) {
                        ind = j;
                        val = tmp(j, i);
                    }
                }

                if (ind != i) {
                    for (index_t j = 0; j != DIM; j++) {
                        val2 = result(i, j);
                        result(i, j) = result(ind, j);
                        result(ind, j) = val2;           /* swap columns */
                        val2 = tmp(i, j);
                        tmp(i, j) = tmp(ind, j);
                        tmp(ind, j) = val2;
                    }
                }

                if (val == 0.0) {
                    return false;
                }

                for (index_t j = 0; j != DIM; j++) {
                    tmp(i, j) /= val;
                    result(i, j) /= val;
                }

                for (index_t j = 0; j != DIM; j++) {
                    if (j == i) {
                        continue;                       /* eliminate column */
                    }
                    val = tmp(j, i);
                    for (index_t k = 0; k != DIM; k++) {
                        tmp(j, k) -= tmp(i, k) * val;
                        result(j, k) -= result(i, k) * val;
                    }
                }
            }

            return true;
        }

        /**
         * \brief Computes the transposed matrix
         * \return the transposed matrix
         */
        matrix_type transpose() const {
            matrix_type result;
            for (index_t i = 0; i < DIM; i++) {
                for (index_t j = 0; j < DIM; j++) {
                    result(i, j) = (*this)(j, i);
                }
            }
            return result;
        }

        /** For interfacing with Fortran, OpenGL etc... */

        /**
         * \brief Gets non-modifiable matrix data
         * \return a const pointer to the first element of the matrix
         */
        inline const FT* data() const {
            return &(coeff_[0][0]);
        }

        /** For interfacing with Fortran, OpenGL etc... */

        /**
         * \brief Gets modifiable matrix data
         * \return a pointer to the first element of the matrix
         */
        inline FT* data() {
            return &(coeff_[0][0]);
        }

        /**
         * \brief Gets the lower triangle of the matrix
         * \details Gets all the coefficients of the matrix under the
         * diagonal (included) to array \p store, Array \p store must be large
         * enough to contain (DIM * (DIM+1))/2 values.
         * \param[in] store an array of at least (DIM * (DIM+1))/2 values
         */
        void get_lower_triangle(FT* store) const {
            for (index_t i = 0; i < DIM; i++) {
                for (index_t j = 0; j <= i; j++) {
                    *store++ = coeff_[i][j];
                }
            }
        }

    private:
        FT coeff_[DIM][DIM];
    };

    /************************************************************************/

    /**
     * \brief Writes a matrix to a stream
     * \details This writes the coefficients of matrix \p m separated by a
     * space character to the output stream \p output.
     * \param[in] output the output stream
     * \param[in] m the matrix to write
     * \return a reference to the output stream \p output
     * \relates Matrix
     */
    template <index_t DIM, class FT>
    inline std::ostream& operator<< (
        std::ostream& output, const Matrix<DIM, FT>& m
        ) {
        const char* sep = "";
        for (index_t i = 0; i < DIM; i++) {
            for (index_t j = 0; j < DIM; j++) {
                output << sep << m(i, j);
                sep = " ";
            }
        }
        return output;
    }

    /**
     * \brief Reads a matrix from a stream
     * \details This reads \p DIM * \p DIM coefficients from the input stream
     * \p input and stores them in matrix \p m
     * \param[in] input the input stream
     * \param[out] m the matrix to read
     * \return a reference to the input stream \p input
     * \relates Matrix
     */
    template <index_t DIM, class FT>
    inline std::istream& operator>> (
        std::istream& input, Matrix<DIM, FT>& m
        ) {
        for (index_t i = 0; i < DIM; i++) {
            for (index_t j = 0; j < DIM; j++) {
                input >> m(i, j);
            }
        }
        return input;
    }

    /************************************************************************/

    /**
     * \brief Multiplies a matrix by a vector
     * \details Multiplies matrix \p M by vector \p x and stores the result in
     * vector y. Vectors \p x and \p y are given as arrays of elements and
     * must at least contain \p DIM elements, otherwise the result is
     * undefined.
     * \param[in] M a \p DIM x \p DIM matrix
     * \param[in] x the input vector
     * \param[in] y the result of the multiplication
     * \tparam FT the type of the matrix elements
     * \tparam DIM the dimension of the matrix
     * \relates Matrix
     */
    template <index_t DIM, class FT> inline
        void mult(const Matrix<DIM, FT>& M, const FT* x, FT* y) {
        for (index_t i = 0; i < DIM; i++) {
            y[i] = 0;
            for (index_t j = 0; j < DIM; j++) {
                y[i] += M(i, j) * x[j];
            }
        }
    }

    /************************************************************************/

    /**
     * \brief Computes a matrix vector product.
     * \param[in] M the matrix
     * \param[in] x the vector
     * \return \p M times \p x
     * \note This function copies the resulting vector, thus it is not
     *  very efficient and should be only used when prototyping.
     */
    template <index_t DIM, class FT> inline
        Vector<DIM, FT> operator*(
            const Matrix<DIM, FT>& M, const Vector<DIM, FT>& x
            ) {
        Vector<DIM, FT> y;
        for (index_t i = 0; i < DIM; i++) {
            y[i] = 0;
            for (index_t j = 0; j < DIM; j++) {
                y[i] += M(i, j) * x[j];
            }
        }
        return y;
    }
    using mat4 = Matrix<4, float>;




    /************************************************************************/

}