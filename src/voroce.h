// APACHE LICENSE, VERSION 2.0
// copyright (c) shinji ogaki

#include "glm/glm.hpp"
#include <array>

namespace voroce
{
	struct Voroce
	{
		// hash function from "Optimized Spatial Hashing for Collision Detection of Deformable Objects"
		// https://matthias-research.github.io/pages/publications/tetraederCollision.pdf

		static const auto PrimeU = 73856093;
		static const auto PrimeV = 19349663;
		static const auto PrimeW = 83492791;

		template <class T>
		static T Hash2D(const T& u, const T& v)
		{
			return (u * PrimeU) ^ (v * PrimeV);
		}

		template <class T>
		static T Hash3D(const T& u, const T& v, const T& w)
		{
			return (u * PrimeU) ^ (v * PrimeV) ^ (w * PrimeW);
		}

		// minstd_rand (TODO: use better hash, do something beter here, fast but a bit ugly...)

		static const auto LCG = 48271;
		template <class T>
		static auto OffsetX(const T& seed)
		{
			return (seed * LCG) / float(0xFFFFFFFF) + 0.5f;
		}

		template <class T>
		static auto OffsetY(const T& seed)
		{
			return (seed * LCG * LCG) / float(0xFFFFFFFF) + 0.5f;
		}

		template <class T>
		static auto OffsetZ(const T& seed)
		{
			return (seed * LCG * LCG * LCG) / float(0xFFFFFFFF) + 0.5f;
		}

		// naive implementation
		static float Evaluate2DRef(const glm::vec2& source)
		{
			// 1 (self) + 20 (neighbours)
			const int us[21] = { 0, 1, 0, -1,  0, 1, -1, -1,  1, 2, 0, -2,  0, 2, 1, -2,  1,  2, -1, -2, -1 };
			const int vs[21] = { 0, 0, 1,  0, -1, 1,  1, -1, -1, 0, 2,  0, -2, 1, 2,  1, -2, -1,  2, -1, -2 };

			const auto flt_x = std::floor(source.x);
			const auto flt_y = std::floor(source.y);
			const auto int_x = (int32_t)flt_x;
			const auto int_y = (int32_t)flt_y;
			const auto origin = glm::vec2(flt_x, flt_y);

			auto worley = std::numeric_limits<float>::max();
			auto result = std::numeric_limits<float>::max();

			for (auto loop = 0; loop < 21; ++loop)
			{
				const auto u = us[loop];
				const auto v = vs[loop];

				const auto hash    = Hash2D(int_x + u, int_y + v);
				const auto offsetX = OffsetX(hash);
				const auto offsetY = OffsetY(hash);
				const auto jitter  = origin + glm::vec2(u + offsetX, v + offsetY);
				const auto tmp     = glm::dot(source - jitter, source - jitter);
				if (worley > tmp)
				{
					worley = tmp;
					result = offsetY;
				}
			}

			//return worley;
			return result;
		}

		// quadrant optimization
		static float Evaluate2DOpt(const glm::vec2& source)
		{
			const auto flt_x = std::floor(source.x);
			const auto flt_y = std::floor(source.y);
			const auto int_x = (int32_t)flt_x;
			const auto int_y = (int32_t)flt_y;
			const auto origin = glm::vec2(flt_x, flt_y);

			// octant selection
			auto quadrant = 0;
			if (0.5f > source.x - origin.x) { quadrant += 1; }
			if (0.5f > source.y - origin.y) { quadrant += 2; }

			const int us[4][16] =
			{
				{0,  1, 0, 1, 0, 1,-1,-1,-1,  2, 2, 0, 1, 2,-1, 2},
				{0, -1,-1, 0,-1, 0, 1, 1, 1, -2,-2,-1, 0,-2, 1,-2},
				{0,  0, 1, 1,-1,-1, 0, 1,-1,  0, 1, 2, 2,-1, 2, 2},
				{0, -1, 0,-1, 1, 1,-1, 0, 1, -1, 0,-2,-2, 1,-2,-2}
			};

			const int vs[4][16] =
			{
				{0,  0, 1, 1,-1,-1, 0, 1,-1,  0, 1, 2, 2,-1, 2, 2},
				{0,  0, 1, 1,-1,-1, 0, 1,-1,  0, 1, 2, 2,-1, 2, 2},
				{0, -1,-1, 0,-1, 0, 1, 1, 1, -2,-2,-1, 0,-2, 1,-2},
				{0, -1,-1, 0,-1, 0, 1, 1, 1, -2,-2,-1, 0,-2, 1,-2}
			};

			// cell id and possible minimum squared distance
			// 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15
			// 0, 0, 0, 0, 1, 1, 1, 1, 2, 4, 4, 4, 4, 5, 5, 8

			auto worley = std::numeric_limits<float>::max();
			auto result = std::numeric_limits<float>::max();

			std::array<float, 5> ranges = { 0, 1, 2,  4,  5 };
			std::array<int  , 5> lowers = { 0, 4, 8,  9, 13 };
			std::array<int  , 5> uppers = { 4, 8, 9, 13, 15 };

			for (auto dist = 0; dist < 5; ++dist)
			{
				if (ranges[dist] < worley * 4)
				{
					for (auto loop = lowers[dist]; loop < uppers[dist]; ++loop)
					{
						const auto u = us[quadrant][loop];
						const auto v = vs[quadrant][loop];
						const auto hash    = Hash2D(int_x + u, int_y + v);
						const auto offsetX = OffsetX(hash);
						const auto offsetY = OffsetY(hash);
						const auto jitter  = origin + glm::vec2(u + offsetX, v + offsetY);
						const auto tmp     = glm::dot(source - jitter, source - jitter);
						if (worley > tmp)
						{
							worley = tmp;
							result = offsetY;
						}
					}
				}
				else
				{
					break;
				}
			}

			//return worley;
			return result;
		}

		// reference implementation
		static float Evaluate3DRef(const glm::vec3& source)
		{
			const auto flt_x = std::floor(source.x);
			const auto flt_y = std::floor(source.y);
			const auto flt_z = std::floor(source.z);
			const auto int_x = (int32_t)flt_x;
			const auto int_y = (int32_t)flt_y;
			const auto int_z = (int32_t)flt_z;
			const auto origin = glm::vec3(flt_x, flt_y, flt_z);

			const auto size = 117;

			const int us[size] = { -1, 0, 1,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-1, 0, 1,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-1, 0, 1,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-1, 0, 1 };
			const int vs[size] = { -2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2,-2,-2,-2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,-2,-2,-2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,-2,-2,-2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,-2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2 };
			const int ws[size] = { -2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 };

			auto worley = std::numeric_limits<float>::max();
			float result;
			for (auto loop = 0; loop < size; ++loop)
			{
				const auto u = us[loop];
				const auto v = vs[loop];
				const auto w = ws[loop];
				const auto hash    = Hash3D(int_x + u, int_y + v, int_z + w);
				const auto offsetX = OffsetX(hash);
				const auto offsetY = OffsetY(hash);
				const auto offsetZ = OffsetZ(hash);
				const auto jitter  = origin + glm::vec3(u + offsetX, v + offsetY, w + offsetZ);
				const auto tmp     = glm::dot(source - jitter, source - jitter);
				if (worley > tmp)
				{
					worley = tmp;
					result = offsetY;
				}
			}

			//return worley;
			return result;
		}

		// octant optimization
		static float Evaluate3DOpt(const glm::vec3& source)
		{
			const auto flt_x = std::floor(source.x);
			const auto flt_y = std::floor(source.y);
			const auto flt_z = std::floor(source.z);
			const auto int_x = (int32_t)flt_x;
			const auto int_y = (int32_t)flt_y;
			const auto int_z = (int32_t)flt_z;
			const auto origin = glm::vec3(flt_x, flt_y, flt_z);

			// octant selection
			auto octant = 0;
			if (0.5f > source.x - origin.x) { octant += 1; }
			if (0.5f > source.y - origin.y) { octant += 2; }
			if (0.5f > source.z - origin.z) { octant += 4; }

			const int us[8][90] =
			{
				{ 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,-1,-1, 0, 1,-1,-1, 0, 1,-1,-1,-1,-1,-1, 2, 2, 0, 1, 2, 2, 0, 1, 0, 1, 0, 1, 2, 2, 0, 1, 2,-1, 2,-1, 0, 1,-1,-1, 2,-1,-1, 2, 2, 2, 2, 0, 1, 0, 1, 0, 1, 2, 0, 1,-2,-2, 0, 1,-2,-2, 2,-1, 0, 1,-1,-1, 0, 1,-2,-2,-1,-2,-1,-2,-1,-1,-2},
				{-1, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1, 0, 1, 1,-1, 0, 1, 1,-1, 0, 1, 1, 1, 1, 1,-2,-2,-1, 0,-2,-2,-1, 0,-1, 0,-1, 0,-2,-2,-1, 0,-2, 1,-2, 1,-1, 0, 1, 1,-2, 1, 1,-2,-2,-2,-2,-1, 0,-1, 0,-1, 0,-2,-1, 0, 2, 2,-1, 0, 2, 2,-2, 1,-1, 0, 1, 1,-1, 0, 2, 2, 1, 2, 1, 2, 1, 1, 2},
				{ 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,-1,-1, 0, 1,-1,-1, 0, 1,-1,-1, 0, 1,-1,-1,-1, 0, 1, 2, 2, 0, 1, 2, 2, 0, 1, 0, 1, 0, 1, 2, 2,-1, 2,-1, 2,-1,-1, 0, 1,-1, 2,-1, 2, 2, 0, 1, 2, 2, 0, 1, 0, 1, 2,-2,-2, 0, 1,-2,-2, 0, 1,-1, 2,-1,-1, 0, 1,-2,-2, 0, 1,-2,-1,-2,-1,-1,-2,-1},
				{-1, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1, 0, 1, 1,-1, 0, 1, 1,-1, 0, 1, 1,-1, 0, 1, 1, 1,-1, 0,-2,-2,-1, 0,-2,-2,-1, 0,-1, 0,-1, 0,-2,-2, 1,-2, 1,-2, 1, 1,-1, 0, 1,-2, 1,-2,-2,-1, 0,-2,-2,-1, 0,-1, 0,-2, 2, 2,-1, 0, 2, 2,-1, 0, 1,-2, 1, 1,-1, 0, 2, 2,-1, 0, 2, 1, 2, 1, 1, 2, 1},
				{ 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,-1,-1, 0, 1,-1,-1, 0, 1, 0, 1,-1,-1, 0, 1,-1,-1,-1, 0, 1, 0, 1, 2, 2, 0, 1, 2, 2, 0, 1, 0, 1,-1,-1, 2,-1, 2,-1, 2, 2, 0, 1,-1, 2,-1, 2, 2, 0, 1, 2, 2, 2,-1, 0, 1,-2,-2, 0, 1,-2,-2, 2, 0, 1, 0, 1,-1,-2,-1,-2, 0, 1,-2,-2, 0, 1,-1,-1,-1,-2,-1},
				{-1, 0,-1, 0,-1, 0,-1, 0,-1, 0, 1, 1,-1, 0, 1, 1,-1, 0,-1, 0, 1, 1,-1, 0, 1, 1, 1,-1, 0,-1, 0,-2,-2,-1, 0,-2,-2,-1, 0,-1, 0, 1, 1,-2, 1,-2, 1,-2,-2,-1, 0, 1,-2, 1,-2,-2,-1, 0,-2,-2,-2, 1,-1, 0, 2, 2,-1, 0, 2, 2,-2,-1, 0,-1, 0, 1, 2, 1, 2,-1, 0, 2, 2,-1, 0, 1, 1, 1, 2, 1},
				{ 0, 1, 0, 1, 0, 1, 0, 1,-1,-1, 0, 1,-1,-1, 0, 1, 0, 1, 0, 1,-1,-1,-1,-1, 0, 1,-1, 0, 1, 0, 1, 0, 1, 2, 2, 0, 1, 2, 2,-1,-1, 0, 1,-1, 2,-1, 2, 0, 1, 2, 2,-1,-1, 2, 0, 1, 2, 2, 2, 2,-1, 2,-2,-2, 0, 1,-2,-2, 0, 1, 2, 0, 1, 0, 1,-2,-1,-2,-1,-2,-2, 0, 1,-1,-1, 0, 1,-2,-1,-1},
				{-1, 0,-1, 0,-1, 0,-1, 0, 1, 1,-1, 0, 1, 1,-1, 0,-1, 0,-1, 0, 1, 1, 1, 1,-1, 0, 1,-1, 0,-1, 0,-1, 0,-2,-2,-1, 0,-2,-2, 1, 1,-1, 0, 1,-2, 1,-2,-1, 0,-2,-2, 1, 1,-2,-1, 0,-2,-2,-2,-2, 1,-2, 2, 2,-1, 0, 2, 2,-1, 0,-2,-1, 0,-1, 0, 2, 1, 2, 1, 2, 2,-1, 0, 1, 1,-1, 0, 2, 1, 1}
			};
			const int vs[8][90] =
			{
				{ 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1,-1,-1, 0, 1,-1,-1, 0, 1,-1,-1, 0, 1,-1,-1,-1, 0, 1, 2, 2, 0, 1, 2, 2, 0, 0, 1, 1, 0, 1, 2, 2,-1, 2,-1, 2,-1,-1, 0, 1,-1, 2,-1, 2, 2, 0, 1, 2, 2, 0, 0, 1, 1, 2,-2,-2, 0, 1,-2,-2, 0, 1,-1, 2,-1,-1, 0, 1,-2,-2, 0, 1,-2,-1,-2,-1,-1,-2,-1},
				{ 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1,-1,-1, 0, 1,-1,-1, 0, 1,-1,-1, 0, 1,-1,-1,-1, 0, 1, 2, 2, 0, 1, 2, 2, 0, 0, 1, 1, 0, 1, 2, 2,-1, 2,-1, 2,-1,-1, 0, 1,-1, 2,-1, 2, 2, 0, 1, 2, 2, 0, 0, 1, 1, 2,-2,-2, 0, 1,-2,-2, 0, 1,-1, 2,-1,-1, 0, 1,-2,-2, 0, 1,-2,-1,-2,-1,-1,-2,-1},
				{-1,-1, 0, 0,-1,-1, 0, 0,-1,-1, 0, 0,-1, 0, 1, 1,-1, 0, 1, 1,-1, 0, 1, 1, 1, 1, 1,-2,-2,-1, 0,-2,-2,-1, 0,-1,-1, 0, 0,-2,-2,-1, 0,-2, 1,-2, 1,-1, 0, 1, 1,-2, 1, 1,-2,-2,-2,-2,-1, 0,-1,-1, 0, 0,-2,-1, 0, 2, 2,-1, 0, 2, 2,-2, 1,-1, 0, 1, 1,-1, 0, 2, 2, 1, 2, 1, 2, 1, 1, 2},
				{-1,-1, 0, 0,-1,-1, 0, 0,-1,-1, 0, 0,-1, 0, 1, 1,-1, 0, 1, 1,-1, 0, 1, 1, 1, 1, 1,-2,-2,-1, 0,-2,-2,-1, 0,-1,-1, 0, 0,-2,-2,-1, 0,-2, 1,-2, 1,-1, 0, 1, 1,-2, 1, 1,-2,-2,-2,-2,-1, 0,-1,-1, 0, 0,-2,-1, 0, 2, 2,-1, 0, 2, 2,-2, 1,-1, 0, 1, 1,-1, 0, 2, 2, 1, 2, 1, 2, 1, 1, 2},
				{ 0, 0, 1, 1, 0, 0, 1, 1,-1,-1, 0, 1,-1,-1, 0, 1, 0, 0, 1, 1,-1,-1,-1,-1, 0, 1,-1, 0, 0, 1, 1, 0, 1, 2, 2, 0, 1, 2, 2,-1,-1, 0, 1,-1, 2,-1, 2, 0, 1, 2, 2,-1,-1, 2, 0, 1, 2, 2, 2, 2,-1, 2,-2,-2, 0, 1,-2,-2, 0, 1, 2, 0, 0, 1, 1,-2,-1,-2,-1,-2,-2, 0, 1,-1,-1, 0, 1,-2,-1,-1},
				{ 0, 0, 1, 1, 0, 0, 1, 1,-1,-1, 0, 1,-1,-1, 0, 1, 0, 0, 1, 1,-1,-1,-1,-1, 0, 1,-1, 0, 0, 1, 1, 0, 1, 2, 2, 0, 1, 2, 2,-1,-1, 0, 1,-1, 2,-1, 2, 0, 1, 2, 2,-1,-1, 2, 0, 1, 2, 2, 2, 2,-1, 2,-2,-2, 0, 1,-2,-2, 0, 1, 2, 0, 0, 1, 1,-2,-1,-2,-1,-2,-2, 0, 1,-1,-1, 0, 1,-2,-1,-1},
				{-1,-1, 0, 0,-1,-1, 0, 0,-1, 0, 1, 1,-1, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1, 0, 1, 1, 1,-1,-1, 0, 0,-2,-2,-1, 0,-2,-2,-1, 0,-1, 0, 1, 1,-2, 1,-2, 1,-2,-2,-1, 0, 1,-2, 1,-2,-2,-1, 0,-2,-2,-2, 1,-1, 0, 2, 2,-1, 0, 2, 2,-2,-1,-1, 0, 0, 1, 2, 1, 2,-1, 0, 2, 2,-1, 0, 1, 1, 1, 2, 1},
				{-1,-1, 0, 0,-1,-1, 0, 0,-1, 0, 1, 1,-1, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1, 0, 1, 1, 1,-1,-1, 0, 0,-2,-2,-1, 0,-2,-2,-1, 0,-1, 0, 1, 1,-2, 1,-2, 1,-2,-2,-1, 0, 1,-2, 1,-2,-2,-1, 0,-2,-2,-2, 1,-1, 0, 2, 2,-1, 0, 2, 2,-2,-1,-1, 0, 0, 1, 2, 1, 2,-1, 0, 2, 2,-1, 0, 1, 1, 1, 2, 1}
			};
			const int ws[8][90] =
			{
				{ 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 0, 1,-1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,-1,-1,-1,-1, 0, 0, 1, 1, 2, 2, 2, 2,-1,-1, 2, 0, 1, 2, 2, 2, 2,-2,-2,-2,-2,-1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2,-2,-2,-2,-2,-1,-1,-1,-1, 0, 0, 1, 1,-2,-1,-1},
				{ 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 0, 1,-1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,-1,-1,-1,-1, 0, 0, 1, 1, 2, 2, 2, 2,-1,-1, 2, 0, 1, 2, 2, 2, 2,-2,-2,-2,-2,-1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2,-2,-2,-2,-2,-1,-1,-1,-1, 0, 0, 1, 1,-2,-1,-1},
				{ 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 0, 1,-1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,-1,-1,-1,-1, 0, 0, 1, 1, 2, 2, 2, 2,-1,-1, 2, 0, 1, 2, 2, 2, 2,-2,-2,-2,-2,-1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2,-2,-2,-2,-2,-1,-1,-1,-1, 0, 0, 1, 1,-2,-1,-1},
				{ 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 0, 1,-1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,-1,-1,-1,-1, 0, 0, 1, 1, 2, 2, 2, 2,-1,-1, 2, 0, 1, 2, 2, 2, 2,-2,-2,-2,-2,-1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2,-2,-2,-2,-2,-1,-1,-1,-1, 0, 0, 1, 1,-2,-1,-1},
				{-1,-1,-1,-1, 0, 0, 0, 0,-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1,-1, 0, 1, 1, 1, 1, 1,-2,-2,-2,-2,-1,-1,-1,-1, 0, 0, 0, 0,-2,-2,-2,-2,-1,-1, 0, 0, 1, 1, 1, 1,-2, 1, 1,-2,-2,-2,-2,-1, 0,-2,-2,-1,-1,-1,-1, 0, 0, 0, 0, 1, 2, 2, 2, 2,-1,-1, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2},
				{-1,-1,-1,-1, 0, 0, 0, 0,-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1,-1, 0, 1, 1, 1, 1, 1,-2,-2,-2,-2,-1,-1,-1,-1, 0, 0, 0, 0,-2,-2,-2,-2,-1,-1, 0, 0, 1, 1, 1, 1,-2, 1, 1,-2,-2,-2,-2,-1, 0,-2,-2,-1,-1,-1,-1, 0, 0, 0, 0, 1, 2, 2, 2, 2,-1,-1, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2},
				{-1,-1,-1,-1, 0, 0, 0, 0,-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1,-1, 0, 1, 1, 1, 1, 1,-2,-2,-2,-2,-1,-1,-1,-1, 0, 0, 0, 0,-2,-2,-2,-2,-1,-1, 0, 0, 1, 1, 1, 1,-2, 1, 1,-2,-2,-2,-2,-1, 0,-2,-2,-1,-1,-1,-1, 0, 0, 0, 0, 1, 2, 2, 2, 2,-1,-1, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2},
				{-1,-1,-1,-1, 0, 0, 0, 0,-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1,-1, 0, 1, 1, 1, 1, 1,-2,-2,-2,-2,-1,-1,-1,-1, 0, 0, 0, 0,-2,-2,-2,-2,-1,-1, 0, 0, 1, 1, 1, 1,-2, 1, 1,-2,-2,-2,-2,-1, 0,-2,-2,-1,-1,-1,-1, 0, 0, 0, 0, 1, 2, 2, 2, 2,-1,-1, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 2}
			};

			// cell id and possible minimum squared distance
			// 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,
			// 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,10,10,10,10,10,10,10,10,10,10,11,11,11,

			auto worley = std::numeric_limits<float>::max();
			auto result = std::numeric_limits<float>::max();

			std::array<float, 11> ranges = { 0,  1,  2,  3,  4,  5,  6,  8,  9, 10, 11 };
			std::array<int,   11> lowers = { 0,  8, 20, 26, 27, 39, 51, 54, 60, 75, 87 };
			std::array<int,   11> uppers = { 8, 20, 26, 27, 39, 51, 54, 60, 75, 87, 90 };

			for (auto dist = 0; dist < 11; ++dist)
			{
				if (ranges[dist] < worley * 4)
				{
					for (auto loop = lowers[dist]; loop < uppers[dist]; ++loop)
					{
						const auto u = us[octant][loop];
						const auto v = vs[octant][loop];
						const auto w = ws[octant][loop];
						const auto hash    = Hash3D(int_x + u, int_y + v, int_z + w);
						const auto offsetX = OffsetX(hash);
						const auto offsetY = OffsetY(hash);
						const auto offsetZ = OffsetZ(hash);
						const auto jitter  = origin + glm::vec3(u + offsetX, v + offsetY, w + offsetZ);
						const auto tmp     = glm::dot(source - jitter, source - jitter);
						if (worley > tmp)
						{
							worley = tmp;
							result = offsetY;
						}
					}
				}
				else
				{
					break;
				}
			}

			//return worley;
			return result;
		}
	};
}