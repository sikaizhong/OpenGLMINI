#pragma once
#include<cassert>
#include<iostream>
#include<fstream>
#include<string>
namespace CGL {
	template <unsigned DIM, class FT>
	class Vector {
	public:
		Vector() {};
		Vector(FT* data) {
			for (unsigned coord = 0; coord < DIM; coord++) {
				data_[coord] = data[coord];
			}
		}

		Vector(FT v0, FT v1, FT v2) {
			data_[0] = v0;
			data_[1] = v1;
			data_[2] = v2;
		}

        Vector(Vector& rh) {
            for (unsigned coord = 0; coord < DIM; coord++) {
                data_[coord] = rh.data_[coord];
            }
        }
        Vector(const Vector& rh) {
            for (unsigned coord = 0; coord < DIM; coord++) {
                data_[coord] = rh.data_[coord];
            }
        }


		//Vector(Vector* rhs)
		const FT& x() const { return data_[0]; }
		const FT& y() const { return data_[1]; }
		const FT& z() const { return data_[2]; }

		
		Vector operator+(const Vector& w) const;
		Vector operator-(const Vector& w) const;

		FT& operator[](unsigned idx) { assert(idx < DIM); return data_[idx]; }
		const FT& operator[](unsigned idx)const { assert(idx < DIM);  return data_[idx]; }
		FT* ptr() { return data_; }
        FT* data() { return data_; }

		unsigned dimension() const {
			return DIM;
		}
		inline FT length2() const {
			FT result = FT(0);
			for (unsigned i = 0; i < DIM; i++) {
				result += data_[i] * data_[i];
			}
			return result;
		}
		inline FT length() const {
			return std::sqrt(length2());
		}

		//friend std::ostream& operator<< (std::ostream& out, const Vector& rhs);
	private:
		FT data_[DIM];
	};

	template <unsigned DIM, class FT>
	inline Vector<DIM, FT>
		Vector<DIM, FT>::operator+(const Vector<DIM, FT>& w) const {
		Vector<DIM, FT> result;
		for (unsigned coord = 0; coord < DIM; coord++) {
			result[coord] = this->data_[coord] + w[coord];
		}
		return result;
	}
	template <unsigned DIM, class FT>
	inline Vector<DIM, FT> Vector<DIM, FT>::operator-(const Vector& w) const {
		Vector<DIM, FT> result;
		for (unsigned coord = 0; coord < DIM; coord++) {
			result[coord] = this->data_[coord]-w[coord];
		}
		return result;
	}



	template <unsigned DIM, class FT>
	inline std::ostream& operator<< (std::ostream& out, const Vector<DIM,FT>& rhs) {
		for (unsigned coord = 0; coord < DIM; coord++) {
			out << rhs[coord] << " ";
		}
		return out;
	}
	template <unsigned DIM, class FT>
	inline Vector<DIM, FT> normalize(const Vector<DIM, FT>& v) {
		FT s = length(v);
		if (s > 1e-30) {
			s = FT(1) / s;
		}
		return s * v;
	}
	template <unsigned DIM, class FT>
	inline FT dot(const Vector<DIM, FT>&v1, const Vector<DIM, FT>&v2){
		FT result = 0;
		for (unsigned i = 0; i < DIM; i++) {
			result += v1[i] * v2[i];
		}
		return result;
	}

	template <unsigned DIM, class FT, class FT2>
	inline Vector<DIM, FT> operator* (FT2 s, const Vector<DIM, FT>& v) {
		Vector<DIM, FT> result;
		for (unsigned i = 0; i < DIM; i++) {
			result[i] = FT(s) * v[i];
		}
		return result;
	}

	template <unsigned DIM, class FT, class FT2>
	inline Vector<DIM, FT> operator* (const Vector<DIM, FT>& v, FT2 s) {
		Vector<DIM, FT> result;
		for (unsigned i = 0; i < DIM; i++) {
			result[i] = FT(s) * v[i];
		}
		return result;
	}

	template <unsigned DIM, class FT>
	inline FT length(const Vector<DIM, FT>& v) {
		return v.length();
	}
    template <class T>
    class Vector<2, T> {
        using index_t = unsigned;
    public:
        static const index_t dim = 2;
        typedef Vector<dim, T> vector_type;
        typedef T value_type;
        Vector() :
            x(0),
            y(0) {
        }
        Vector(T x_in, T y_in) :
            x(x_in),
            y(y_in) {
        }

        template <class T2>
        explicit Vector(const Vector<dim, T2>& v) :
            x(v.x),
            y(v.y) {
        }

        template <class T2>
        explicit Vector(const T2* v) :
            x(v[0]),
            y(v[1]) {
        }

        inline T length2() const {
            return x * x + y * y;
        }

        inline T length() const {
            return sqrt(x * x + y * y);
        }

        inline T distance2(const vector_type& v) const {
            T dx = v.x - x;
            T dy = v.y - y;
            return dx * dx + dy * dy;
        }

        inline T distance(const vector_type& v) const {
            return sqrt(distance2(v));
        }

        inline vector_type& operator+= (const vector_type& v) {
            x += v.x;
            y += v.y;
            return *this;
        }

        inline vector_type& operator-= (const vector_type& v) {
            x -= v.x;
            y -= v.y;
            return *this;
        }

        template <class T2>
        inline vector_type& operator*= (T2 s) {
            x *= T(s);
            y *= T(s);
            return *this;
        }

        template <class T2>
        inline vector_type& operator/= (T2 s) {
            x /= T(s);
            y /= T(s);
            return *this;
        }

        inline vector_type operator+ (const vector_type& v) const {
            return vector_type(x + v.x, y + v.y);
        }

        inline vector_type operator- (const vector_type& v) const {
            return vector_type(x - v.x, y - v.y);
        }

        template <class T2>
        inline vector_type operator* (T2 s) const {
            return vector_type(x * T(s), y * T(s));
        }

        template <class T2>
        inline vector_type operator/ (T2 s) const {
            return vector_type(x / T(s), y / T(s));
        }

        inline vector_type operator- () const {
            return vector_type(-x, -y);
        }

        index_t dimension() const {
            return dim;
        }

        T* data() {
            return &x;
        }
        const T* data() const {
            return &x;
        }

        inline T& operator[] (index_t i) {
            return data()[i];
        }

        inline const T& operator[] (index_t i) const {
            return data()[i];
        }
        T x, y;
    };

    template <class T>
    inline T dot(
        const Vector<2, T>& v1, const Vector<2, T>& v2
    ) {
        return v1.x * v2.x + v1.y * v2.y;
    }
    template <class T>
    inline T det(
        const Vector<2, T>& v1, const Vector<2, T>& v2
    ) {
        return v1.x * v2.y - v1.y * v2.x;
    }
    template <class T2, class T>
    inline Vector<2, T> operator* (
        T2 s, const Vector<2, T>& v
        ) {
        return vecng<2, T>(T(s) * v.x, T(s) * v.y);
    }
    template <class T>
    class Vector<3,T> {
    public:
        static const unsigned dim = 3;
        using  vector_type= Vector<3, T> ;
        using value_type = T;
        Vector() : x(0),y(0), z(0) { }
        Vector(T x_in, T y_in, T z_in) :x(x_in),y(y_in),z(z_in) {}
        template <class T2>
        explicit Vector(const Vector<dim, T2>& v) :x(v.x), y(v.y),z(v.z) {}
        template <class T2>
        explicit Vector(const T2* v) :x(v[0]),y(v[1]),z(v[2]) {}
        inline T length2() const {return x * x + y * y + z * z;}
        inline T length() const {return sqrt(x * x + y * y + z * z);}
        inline T distance2(const vector_type& v) const {
            T dx = v.x - x;
            T dy = v.y - y;
            T dz = v.z - z;
            return dx * dx + dy * dy + dz * dz;
        }
        inline T distance(const vector_type& v) const {
            return sqrt(distance2(v));
        }
        inline vector_type& operator+= (const vector_type& v) {
            x += v.x;
            y += v.y;
            z += v.z;
            return *this;
        }
        inline vector_type& operator-= (const vector_type& v) {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            return *this;
        }
        template <class T2>
        inline vector_type& operator*= (T2 s) {
            x *= T(s);
            y *= T(s);
            z *= T(s);
            return *this;
        }
        template <class T2>
        inline vector_type& operator/= (T2 s) {
            x /= T(s);
            y /= T(s);
            z /= T(s);
            return *this;
        }
        inline vector_type operator+ (const vector_type& v) const {
            return vector_type(x + v.x, y + v.y, z + v.z);
        }
        inline vector_type operator- (const vector_type& v) const {
            return vector_type(x - v.x, y - v.y, z - v.z);
        }
        template <class T2>
        inline vector_type operator* (T2 s) const {
            return vector_type(x * T(s), y * T(s), z * T(s));
        }
        template <class T2>
        inline vector_type operator/ (T2 s) const {
            return vector_type(x / T(s), y / T(s), z / T(s));
        }
        inline vector_type operator- () const {
            return vector_type(-x, -y, -z);
        }
        unsigned dimension() const {
            return dim;
        }
        T* data() {return &x;}

        const T* data() const {return &x;}
        inline T& operator[] (unsigned i) {
            return data()[i];
        }
        inline const T& operator[] (unsigned i) const {
            return data()[i];
        }
        T x, y, z;
    };
    template <class T>
    inline T dot(
        const Vector<3,T>& v1, const Vector<3,T>& v2
    ) {
        return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
    }
    template <class T>
    inline Vector<3,T> cross(
        const Vector<3,T>& v1, const Vector<3,T>& v2
    ) {
        return Vector<3,T>(
            v1.y * v2.z - v1.z * v2.y,
            v1.z * v2.x - v1.x * v2.z,
            v1.x * v2.y - v1.y * v2.x
        );
    }

    template <class T2, class T>
    inline Vector<3,T> operator* (
        T2 s, const Vector<3,T>& v
        ) {
        return Vector<3,T>(T(s) * v.x, T(s) * v.y, T(s) * v.z);
    }



    template <class T>
    class Vector<4, T> {
    public:
        static const unsigned dim = 4;
        using  vector_type = Vector<4, T>;
        using value_type = T;
        Vector() : x(0), y(0), z(0),w(0) { }
        Vector(T x_in, T y_in, T z_in, T w_in) :x(x_in), y(y_in), z(z_in),w(w_in) {}
        template <class T2>
        explicit Vector(const Vector<dim, T2>& v) :x(v.x), y(v.y), z(v.z),w(v.w) {}
        template <class T2>
        explicit Vector(const T2* v) :x(v[0]), y(v[1]), z(v[2]),w(v[3]) {}
        inline T length2() const { return x * x + y * y + z * z+w*w; }
        inline T length() const { return sqrt(x * x + y * y + z * z+w*w); }
        inline T distance2(const vector_type& v) const {
            T dx = v.x - x;
            T dy = v.y - y;
            T dz = v.z - z;
            T dw = v.w - w;
            return dx * dx + dy * dy + dz * dz+dw*dw;
        }
        inline T distance(const vector_type& v) const {
            return sqrt(distance2(v));
        }
        inline vector_type& operator+= (const vector_type& v) {
            x += v.x;
            y += v.y;
            z += v.z;
            w += v.w;
            return *this;
        }
        inline vector_type& operator-= (const vector_type& v) {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            w -= v.w;
            return *this;
        }
        template <class T2>
        inline vector_type& operator*= (T2 s) {
            x *= T(s);
            y *= T(s);
            z *= T(s);
            w *= T(s);
            return *this;
        }
        template <class T2>
        inline vector_type& operator/= (T2 s) {
            x /= T(s);
            y /= T(s);
            z /= T(s);
            w /= T(s);
            return *this;
        }
        inline vector_type operator+ (const vector_type& v) const {
            return vector_type(x + v.x, y + v.y, z + v.z, w+v.w);
        }
        inline vector_type operator- (const vector_type& v) const {
            return vector_type(x - v.x, y - v.y, z - v.z, w-v.w);
        }
        template <class T2>
        inline vector_type operator* (T2 s) const {
            return vector_type(x * T(s), y * T(s), z * T(s), w*T(s));
        }
        template <class T2>
        inline vector_type operator/ (T2 s) const {
            return vector_type(x / T(s), y / T(s), z / T(s),w/T(s));
        }
       inline vector_type operator- () const {
            return vector_type(-x, -y, -z,-w);
        }
        unsigned dimension() const {
            return dim;
        }
        T* data() { return &x; }

        const T* data() const { return &x; }
        inline T& operator[] (unsigned i) {
            return data()[i];
        }
        inline const T& operator[] (unsigned i) const {
            return data()[i];
        }
        T x, y, z,w;
    };










    using vec2 = Vector<2, float>;
    using ivec2 = Vector<2, int>;

    using vec3 = Vector<3, float>;
    using vec4 = Vector<4, float>;
}