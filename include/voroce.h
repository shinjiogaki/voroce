// APACHE LICENSE, VERSION 2.0
// copyright (c) shinji ogaki

#include "glm/glm.hpp"
#include <array>
#include <functional>

struct Voroce
{
	// hash function from "Optimized Spatial Hashing for Collision Detection of Deformable Objects"
	// https://matthias-research.github.io/pages/publications/tetraederCollision.pdf

	static const auto PrimeU = 73856093;
	static const auto PrimeV = 19349663;
	static const auto PrimeW = 83492791;
	static const auto LCG    = 48271;

	static const int32_t Hash2DLowQuality(const int32_t u, const int32_t v)
	{
		return ((u * PrimeU) ^ (v * PrimeV)) * LCG;
	}

	static const int32_t Hash3DLowQuality(const int32_t u, const int32_t v, const int32_t w)
	{
		return ((u * PrimeU) ^ (v * PrimeV) ^ (w * PrimeW)) * LCG;
	}

	// minstd_rand (TODO: use better hash, do something beter here, fast but a bit ugly...)

	static auto OffsetX(const int32_t seed, const float jitter)
	{
		return ((seed) / float(0xFFFFFFFF)) * jitter + 0.5f;
	}

	static auto OffsetY(const int32_t seed, const float jitter)
	{
		return ((seed * LCG) / float(0xFFFFFFFF)) * jitter + 0.5f;
	}

	static auto OffsetZ(const int32_t seed, const float jitter)
	{
		return ((seed * LCG * LCG) / float(0xFFFFFFFF)) * jitter + 0.5f;
	}

	// naive implementation
	static std::pair<int32_t, float> Evaluate2DRef(const glm::vec2& source, const int32_t (*hash_2d)(const int32_t, const int32_t), const float jitter = 1.0f)
	{
		assert(0.0f <= jitter && jitter <= 1.0f);

		const auto flt_x = std::floor(source.x);
		const auto flt_y = std::floor(source.y);
		const auto int_x = int32_t(flt_x);
		const auto int_y = int32_t(flt_y);
		const auto origin = glm::vec2(flt_x, flt_y);

		auto sq_dist = std::numeric_limits<float>::max();
		auto cell_id = 0;

		// 1 (self) + 20 (neighbours)
		const std::array<int32_t, 21> us = { 0, 1, 0, -1,  0, 1, -1, -1,  1, 2, 0, -2,  0, 2, 1, -2,  1,  2, -1, -2, -1 };
		const std::array<int32_t, 21> vs = { 0, 0, 1,  0, -1, 1,  1, -1, -1, 0, 2,  0, -2, 1, 2,  1, -2, -1,  2, -1, -2 };
		for (auto loop = 0; loop < 21; ++loop)
		{
			const auto u = us[loop];
			const auto v = vs[loop];

			const auto hash    = hash_2d(int_x + u, int_y + v);
			const auto offsetX = OffsetX(hash, jitter);
			const auto offsetY = OffsetY(hash, jitter);
			const auto sample  = origin + glm::vec2(u + offsetX, v + offsetY);
			const auto tmp     = glm::dot(source - sample, source - sample);
			if (sq_dist > tmp)
			{
				sq_dist = tmp;
				cell_id = hash;
			}
		}

		return std::make_pair(cell_id, sq_dist);
	}

	// quadrant optimization
	static std::pair<int, float> Evaluate2DOpt(const glm::vec2& source, const int32_t (*hash_2d)(const int32_t, const int32_t), const float jitter = 1.0f)
	{
		assert(0.0f <= jitter && jitter <= 1.0f);

		const auto flt_x = std::floor(source.x);
		const auto flt_y = std::floor(source.y);
		const auto int_x = int32_t(flt_x);
		const auto int_y = int32_t(flt_y);
		const auto origin = glm::vec2(flt_x, flt_y);

		// early termination
		if (0.0f == jitter)
		{
			const auto hash   = hash_2d(int_x, int_y);
			const auto sample = origin + 0.5f;
			return std::make_pair(hash, glm::dot(source - sample, source - sample));
		}

		// octant selection
		auto quadrant = 0;
		if (0.5f > source.x - origin.x) { quadrant += 1; }
		if (0.5f > source.y - origin.y) { quadrant += 2; }

		const int32_t us[4][15] =
		{
			{0,  1, 0, 1, 0, 1,-1,-1,-1,  2, 2, 0, 1, 2,-1},
			{0, -1,-1, 0,-1, 0, 1, 1, 1, -2,-2,-1, 0,-2, 1},
			{0,  0, 1, 1,-1,-1, 0, 1,-1,  0, 1, 2, 2,-1, 2},
			{0, -1, 0,-1, 1, 1,-1, 0, 1, -1, 0,-2,-2, 1,-2}
		};

		const int32_t vs[4][15] =
		{
			{0,  0, 1, 1,-1,-1, 0, 1,-1,  0, 1, 2, 2,-1, 2},
			{0,  0, 1, 1,-1,-1, 0, 1,-1,  0, 1, 2, 2,-1, 2},
			{0, -1,-1, 0,-1, 0, 1, 1, 1, -2,-2,-1, 0,-2, 1},
			{0, -1,-1, 0,-1, 0, 1, 1, 1, -2,-2,-1, 0,-2, 1}
		};

		// cell id and possible minimum squared distance
		// 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14
		// 0, 0, 0, 0, 1, 1, 1, 1, 2, 4, 4, 4, 4, 5, 5

		auto sq_dist = std::numeric_limits<float>::max();
		auto cell_id = 0;

		std::array<  float, 5> ranges = { 0, 1, 2,  4,  5 };
		std::array<int32_t, 5> lowers = { 0, 4, 8,  9, 13 };
		std::array<int32_t, 5> uppers = { 4, 8, 9, 13, 15 };

		for (auto dist = 0; dist < 5; ++dist)
		{
			if (ranges[dist] < sq_dist * 4)
			{
				for (auto loop = lowers[dist]; loop < uppers[dist]; ++loop)
				{
					const auto u = us[quadrant][loop];
					const auto v = vs[quadrant][loop];
					const auto hash    = hash_2d(int_x + u, int_y + v);
					const auto offsetX = OffsetX(hash, jitter);
					const auto offsetY = OffsetY(hash, jitter);
					const auto sample  = origin + glm::vec2(u + offsetX, v + offsetY);
					const auto tmp     = glm::dot(source - sample, source - sample);
					if (sq_dist > tmp)
					{
						sq_dist = tmp;
						cell_id = hash;
					}
				}
			}
			else
			{
				break;
			}
		}

		return std::make_pair(cell_id, sq_dist);
	}

	// reference implementation
	static std::pair<int, float> Evaluate3DRef(const glm::vec3& source, const int32_t (*hash_3d)(const int32_t, const int32_t, const int32_t), const float jitter = 1.0f)
	{
		assert(0.0f <= jitter && jitter <= 1.0f);

		const auto flt_x = std::floor(source.x);
		const auto flt_y = std::floor(source.y);
		const auto flt_z = std::floor(source.z);
		const auto int_x = int32_t(flt_x);
		const auto int_y = int32_t(flt_y);
		const auto int_z = int32_t(flt_z);
		const auto origin = glm::vec3(flt_x, flt_y, flt_z);

		const auto size = 117;

		const std::array<int32_t, size> us = { -1, 0, 1,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-1, 0, 1,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-1, 0, 1,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-1, 0, 1 };
		const std::array<int32_t, size> vs = { -2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2,-2,-2,-2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,-2,-2,-2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,-2,-2,-2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,-2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2 };
		const std::array<int32_t, size> ws = { -2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 };

		auto sq_dist = std::numeric_limits<float>::max();
		auto cell_id = 0;
		for (auto loop = 0; loop < size; ++loop)
		{
			const auto u = us[loop];
			const auto v = vs[loop];
			const auto w = ws[loop];
			const auto hash    = hash_3d(int_x + u, int_y + v, int_z + w);
			const auto offsetX = OffsetX(hash, jitter);
			const auto offsetY = OffsetY(hash, jitter);
			const auto offsetZ = OffsetZ(hash, jitter);
			const auto sample  = origin + glm::vec3(u + offsetX, v + offsetY, w + offsetZ);
			const auto tmp     = glm::dot(source - sample, source - sample);
			if (sq_dist > tmp)
			{
				sq_dist = tmp;
				cell_id = hash;
			}
		}

		return std::make_pair(cell_id, sq_dist);
	}

	// octant optimization
	static std::pair<int, float> Evaluate3DOpt(const glm::vec3& source, const int32_t (*hash_3d)(const int32_t, const int32_t, const int32_t), const float jitter = 1.0f)
	{
		assert(0.0f <= jitter && jitter <= 1.0f);

		const auto flt_x = std::floor(source.x);
		const auto flt_y = std::floor(source.y);
		const auto flt_z = std::floor(source.z);
		const auto int_x = int32_t(flt_x);
		const auto int_y = int32_t(flt_y);
		const auto int_z = int32_t(flt_z);
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

		auto sq_dist = std::numeric_limits<float>::max();
		auto cell_id = 0;
		std::array<  float, 11> ranges = { 0,  1,  2,  3,  4,  5,  6,  8,  9, 10, 11 };
		std::array<int32_t, 11> lowers = { 0,  8, 20, 26, 27, 39, 51, 54, 60, 75, 87 };
		std::array<int32_t, 11> uppers = { 8, 20, 26, 27, 39, 51, 54, 60, 75, 87, 90 };

		for (auto dist = 0; dist < 11; ++dist)
		{
			if (ranges[dist] < sq_dist * 4)
			{
				for (auto loop = lowers[dist]; loop < uppers[dist]; ++loop)
				{
					const auto u = us[octant][loop];
					const auto v = vs[octant][loop];
					const auto w = ws[octant][loop];
					const auto hash    = hash_3d(int_x + u, int_y + v, int_z + w);
					const auto offsetX = OffsetX(hash, jitter);
					const auto offsetY = OffsetY(hash, jitter);
					const auto offsetZ = OffsetZ(hash, jitter);
					const auto sample  = origin + glm::vec3(u + offsetX, v + offsetY, w + offsetZ);
					const auto tmp     = glm::dot(source - sample, source - sample);
					if (sq_dist > tmp)
					{
						sq_dist = tmp;
						cell_id = hash;
					}
				}
			}
			else
			{
				break;
			}
		}

		return std::make_pair(cell_id, sq_dist);
	}
};
