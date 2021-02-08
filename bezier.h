#ifndef BEZIER_H
#define BEZIER_H
#include <iostream>
#include <functional>
#include <cmath>
#include <vector>
#include <memory>
namespace bezier {
	
	namespace types {
		using real_t = double;
		
		using node_index_t = unsigned int;
		
		class point_2d {
		public:
			const real_t X, Y;
			
			point_2d (const real_t& c_x, const real_t& c_y) :
					X(c_x), Y(c_y) {}
			
			friend point_2d operator+(const point_2d&, const point_2d&);
			
			friend point_2d operator*(const point_2d&, const real_t&);
			
			friend std::ostream& operator<<(std::ostream&, const point_2d&);
		};
		point_2d operator+(const point_2d& p_one, const point_2d& p_two) {
		return point_2d(p_one.X + p_two.X, 
							p_one.Y + p_two.Y);
		}

		point_2d operator*(const point_2d& p, const real_t& scalar) {
			return point_2d(p.X * scalar, p.Y * scalar);
		}

		point_2d operator*(const real_t& scalar, const point_2d& p) {
			return p*scalar;
		}

		std::ostream& operator<<(std::ostream& os, const point_2d& p) {
			os << "(" << p.X << ", " << p.Y << ")";
			return os;
		}
	}
	
	namespace constants {
		const int NUM_OF_CUBIC_BEZIER_NODES = 4;
		
		const double ARC = 4.*(sqrt(2.) - 1.)/3.;
	}
	
	using namespace types;
	using namespace constants;
	
	using curve_t = typename std::function<point_2d(const node_index_t&)>;
	
	template<typename T>
	std::function<T(const node_index_t&)> func(T x) {
		return [x](const node_index_t& index) {
			if(index == 0) {
				return x;
			}
			else {
				throw std::out_of_range("a curve node index is out of range");
			}
		};
	}
	
	template<typename T, typename... Args>
	std::function<T(const node_index_t&)> func(T x, Args... args) {
		return [x, args...](const node_index_t& index) {
			if(index == 0) {
				return x;
			}
			else {
				return func(args...)(index - 1);
			}
		};
	}
	
	auto Cup() {
		point_2d points[NUM_OF_CUBIC_BEZIER_NODES]{{-1., 1.}, {-1., -1.}, {1., -1.}, {1., 1.}};
		return func(points[0], points[1], points[2], points[3]);
	}
	
	auto Cap() {
		point_2d points[NUM_OF_CUBIC_BEZIER_NODES]{{-1., -1.}, {-1., 1.}, {1., 1.}, {1., -1.}};
		return func(points[0], points[1], points[2], points[3]);
	}
	
	auto ConvexArc() {
		point_2d points[NUM_OF_CUBIC_BEZIER_NODES]{{0., 1.}, {ARC, 1.}, {1., ARC}, {1., 0}};
		return func(points[0], points[1], points[2], points[3]);
	}
	
	auto ConcaveArc() {
		point_2d points[NUM_OF_CUBIC_BEZIER_NODES]{{0., 1.}, {0, 1. - ARC}, {1. - ARC, 0.}, {1., 0}};
		return func(points[0], points[1], points[2], points[3]);
	}
	
	auto LineSegment(const point_2d& p, const point_2d& q) {
		return func(p, p, q, q);
	}
	
	curve_t Concatenate(curve_t f) {
		return f;
	}
	
	curve_t Concatenate(curve_t f_1, curve_t f_2) {
		return [f_1, f_2](const node_index_t& index) {
			if(index < NUM_OF_CUBIC_BEZIER_NODES) {
				return f_1(index);
			}
			else {
				return f_2(index - NUM_OF_CUBIC_BEZIER_NODES);
			}
		};
	}
	
	template<typename... Args>
	curve_t Concatenate(curve_t f_1, Args... args) {
		return Concatenate(f_1, Concatenate(args...));
	}
	
	auto MovePoint(curve_t f, const node_index_t& mod_index, real_t x, real_t y) {
		return [f, mod_index, x, y](const node_index_t& index) {
			if(index == mod_index) {
				return point_2d(x, y) + f(index);
			}
			else {
				return f(index);
			}
		};
	}
	
	auto Rotate(curve_t f, const double& angle) {
			return [f, angle](const node_index_t& index) {
				point_2d ret_val = f(index);
				double radians = angle*M_PI/180;
				double sine = sin(radians);
				double cosine = cos(radians);
				real_t x = ret_val.X;
				real_t y = ret_val.Y;
				return point_2d(x * cosine - y * sine, 
								x * sine + y * cosine);
			};
	}
	
	auto Scale(curve_t f, const double& x_scale, const double& y_scale) {
		return [f, x_scale, y_scale](const node_index_t& index) {
			point_2d ret_val = f(index);
			real_t x = ret_val.X;
			real_t y = ret_val.Y;
			return point_2d(x * x_scale, y * y_scale);
		};
	}
	
	auto Translate(curve_t f, const double& x_move, const double& y_move) {
		return [f, x_move, y_move](const node_index_t& index) {
			return f(index) + point_2d(x_move, y_move);
		};
	}
	
	class P3CurvePlotter {
	private:
		std::vector<std::vector<bool>> bitmap;
		double step;
		real_t abs(const real_t& real) {
			if(real < 0.) {
				 return -real;
			}
			else {
				return real;
			}
		}
		
		bool is_in_printed_box(const point_2d& p) {
			return abs(p.X) <= 1. && abs(p.Y) <= 1.;
		}
		
		void deCasteljau(const std::vector<point_2d>& control_points,
										  std::vector<point_2d>& points_on_curve) {
			for(double t = 0; t <= 1; t += step) {
				std::vector<point_2d> b(control_points);
				int ptr = 0;
				for(int i = 0; i < NUM_OF_CUBIC_BEZIER_NODES; ++i) {
					for(int j = 0; j < NUM_OF_CUBIC_BEZIER_NODES - i - 1; ++j) {
						b.push_back((1. - t) * b[ptr] + t * b[ptr + 1]);
						++ptr; 
					}
					++ptr;
				}
				points_on_curve.push_back(b[b.size() - 1]);
			}
		}
		
		std::vector<point_2d> calculate_points(curve_t f, const unsigned int& num_of_segments) {
			std::vector<point_2d> points_on_curve;
			for(unsigned int i = 0; i < num_of_segments; ++i) {
				std::vector<point_2d> control_points;
				for(unsigned int j = 0; j < NUM_OF_CUBIC_BEZIER_NODES; ++j) {
					control_points.push_back(f(i * NUM_OF_CUBIC_BEZIER_NODES + j));
				}
				deCasteljau(control_points, points_on_curve);
			}
			return points_on_curve;
		}
		
	public:
		P3CurvePlotter(curve_t f, 
					   const unsigned int& num_of_segments = 1, 
					   const unsigned int& resolution = 80) {
			bitmap = std::vector<std::vector<bool>>
					(resolution, std::vector<bool>(resolution, 0));
			
			step = pow(2./((double)resolution), 2);
			
			auto points_on_curve =
						calculate_points(f, num_of_segments);
			
			for(auto point : points_on_curve) {
				if(is_in_printed_box(point)) {
					real_t x = point.X + 1.; //[0, 2]
					real_t y = 1. - point.Y; //[0, 2]
					x *= (double)resolution/2;
					y *= (double)resolution/2;
					unsigned int normalized_x = x;
					unsigned int normalized_y = y;
					if(normalized_x == resolution) --normalized_x;
					if(normalized_y == resolution) --normalized_y;
					bitmap[normalized_y][normalized_x] = 1;
				}
			}
		}
		
		void Print(std::ostream& os = std::cout, const char curve_sign = '*', const char background_sign = ' ') const {
			for(auto row : bitmap) {
				for(auto bit : row) {
					if(bit) {
						os << curve_sign;
					}
					else {
						os << background_sign;
					}
				}
				os << std::endl;
			}
		}
		
		types::point_2d operator() (curve_t f, double t, int segment_num) {
			int first_id = 4*(segment_num-1);
			std::vector<point_2d> p({f(first_id), f(first_id + 1), f(first_id + 2), f(first_id + 3)});
			std::vector<point_2d> returned;
			deCasteljau(returned, p);
			return returned[0];
		} 
		 
	};
}
#endif
