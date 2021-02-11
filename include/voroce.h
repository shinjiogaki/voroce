// APACHE LICENSE, VERSION 2.0
// copyright (c) shinji ogaki

#include "glm/glm.hpp"
#include <array>

struct Voroce
{
	// hash function from "Optimized Spatial Hashing for Collision Detection of Deformable Objects"
	// https://matthias-research.github.io/pages/publications/tetraederCollision.pdf

	static const auto PrimeU = 73856093;
	static const auto PrimeV = 13467053;
	static const auto PrimeW = 83492791;
	static const auto PrimeT = 23469181;
	static const auto LCG    = 48271;

	static const int32_t Hash2DLowQuality(const int32_t u, const int32_t v)
	{
		return ((u * PrimeU) ^ (v * PrimeV)) * LCG;
	}

	static const int32_t Hash3DLowQuality(const int32_t u, const int32_t v, const int32_t w)
	{
		return ((u * PrimeU) ^ (v * PrimeV) ^ (w * PrimeW)) * LCG;
	}

	static const int32_t Hash4DLowQuality(const int32_t u, const int32_t v, const int32_t w, const int32_t t)
	{
		return ((u * PrimeU) ^ (v * PrimeV) ^ (w * PrimeW) ^ (t * PrimeT)) * LCG;
	}

	// minstd_rand (TODO: use better hash, do something beter here, fast but a bit ugly...)

	static auto OffsetX(const int32_t seed)
	{
		return ((seed) / float(0xFFFFFFFF)) + 0.5f;
	}

	static auto OffsetY(const int32_t seed)
	{
		return ((seed * LCG) / float(0xFFFFFFFF)) + 0.5f;
	}

	static auto OffsetZ(const int32_t seed)
	{
		return (((seed * LCG) * LCG) / float(0xFFFFFFFF)) + 0.5f;
	}

	static auto OffsetT(const int32_t seed)
	{
		return ((((seed * LCG) * LCG) * LCG) / float(0xFFFFFFFF)) + 0.5f;
	}

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
		return (((seed * LCG) * LCG) / float(0xFFFFFFFF)) * jitter + 0.5f;
	}

	static auto OffsetT(const int32_t seed, const float jitter)
	{
		return ((((seed * LCG) * LCG) * LCG) / float(0xFFFFFFFF)) * jitter + 0.5f;
	}

	// naive implementation
	static std::pair<int32_t, float> Evaluate2DRef(const glm::vec2& source, const int32_t (*my_hash)(const int32_t, const int32_t), const float jitter = 1.0f)
	{
		assert(0.0f <= jitter && jitter <= 1.0f);

		const auto flt_x  = std::floor(source.x);
		const auto flt_y  = std::floor(source.y);
		const auto int_x  = int32_t(flt_x);
		const auto int_y  = int32_t(flt_y);
		const auto origin = glm::vec2(flt_x, flt_y);

		auto sq_dist = std::numeric_limits<float>::max();
		auto cell_id = 0;

		// 1 (self) + 20 (neighbours)
		const auto size = 21;
		const std::array<int32_t, size> us = { 0, 0,-1, 1, 0,-1, 1,-1, 1, 0,-2, 2, 0,-1, 1,-2, 2,-2, 2,-1, 1 };
		const std::array<int32_t, size> vs = { 0,-1, 0, 0, 1,-1,-1, 1, 1,-2, 0, 0, 2,-2,-2,-1,-1, 1, 1, 2, 2 };

		for (auto loop = 0; loop < size; ++loop)
		{
			const auto u = us[loop];
			const auto v = vs[loop];

			const auto hash    = my_hash(int_x + u, int_y + v);
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
	static std::pair<int, float> Evaluate2DOpt(const glm::vec2& source, const int32_t (*my_hash)(const int32_t, const int32_t), const float jitter = 1.0f)
	{
		assert(0.0f <= jitter && jitter <= 1.0f);

		const auto flt_x  = std::floor(source.x);
		const auto flt_y  = std::floor(source.y);
		const auto int_x  = int32_t(flt_x);
		const auto int_y  = int32_t(flt_y);
		const auto origin = glm::vec2(flt_x, flt_y);

		// early termination
		if (0.0f == jitter)
		{
			const auto hash   = my_hash(int_x, int_y);
			const auto sample = origin + 0.5f;
			return std::make_pair(hash, glm::dot(source - sample, source - sample));
		}

		// quadrant selection
		auto quadrant = 0;
		if (0.5f > source.x - origin.x) { quadrant += 1; }
		if (0.5f > source.y - origin.y) { quadrant += 2; }

		const int32_t us[4][13] =
		{
			{ 0, 1, 0, 1, 0,-1, 1,-1,-1, 2, 0, 2, 1},
			{ 0,-1, 0,-1, 0, 1,-1, 1, 1,-2, 0,-2,-1},
			{ 0, 0, 1, 1,-1, 0,-1, 1,-1, 0, 2, 1, 2},
			{ 0, 0,-1,-1, 1, 0, 1,-1, 1, 0,-2,-1,-2},
		};

		const int32_t vs[4][13] =
		{
			{ 0, 0, 1, 1,-1, 0,-1, 1,-1, 0, 2, 1, 2},
			{ 0, 0, 1, 1,-1, 0,-1, 1,-1, 0, 2, 1, 2},
			{ 0,-1, 0,-1, 0, 1,-1, 1, 1,-2, 0,-2,-1},
			{ 0,-1, 0,-1, 0, 1,-1, 1, 1,-2, 0,-2,-1},
		};

		auto sq_dist = std::numeric_limits<float>::max();
		auto cell_id = 0;

		if (0.5f >= jitter)
		{
			for (auto loop = 0; loop < 4; ++loop)
			{
				const auto u       = us[quadrant][loop];
				const auto v       = vs[quadrant][loop];
				const auto hash    = my_hash(int_x + u, int_y + v);
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

		// cell id and possible minimum squared distance
		// 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12
		// 0, 0, 0, 0, 1, 1, 1, 1, 2, 4, 4, 4, 4
		const std::array<int32_t, 5> slices = { 0, 4, 8, 9, 13 };
		const std::array<int32_t, 4> ranges = { 0, 1, 2, 4 };

		for (auto dist = 0; dist < 4; ++dist)
		{
			if (ranges[dist] < sq_dist * 4)
			{
				for (auto loop = slices[dist]; loop < slices[dist + 1]; ++loop)
				{
					const auto u       = us[quadrant][loop];
					const auto v       = vs[quadrant][loop];
					const auto hash    = my_hash(int_x + u, int_y + v);
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

	// quadrant optimization
	int32_t   cache2d_cell_id  = 0xFFFFFFFF;
	int32_t   cache2d_quadrant = 0;
	int32_t   cache2d_counter  = 0;
	int32_t   cache2d_cell_ids[16];
	glm::vec2 cache2d_samples [16];
	std::pair<int, float> Evaluate2DCache(const glm::vec2& source, const int32_t(*my_hash)(const int32_t, const int32_t), const float jitter = 1.0f)
	{
		assert(0.0f <= jitter && jitter <= 1.0f);

		const auto flt_x  = std::floor(source.x);
		const auto flt_y  = std::floor(source.y);
		const auto int_x  = int32_t(flt_x);
		const auto int_y  = int32_t(flt_y);
		const auto origin = glm::vec2(flt_x, flt_y);

		// early termination
		if (0.0f == jitter)
		{
			const auto hash   = my_hash(int_x, int_y);
			const auto sample = origin + 0.5f;
			return std::make_pair(hash, glm::dot(source - sample, source - sample));
		}

		// quadrant selection
		auto quadrant = 0;
		if (0.5f > source.x - origin.x) { quadrant += 1; }
		if (0.5f > source.y - origin.y) { quadrant += 2; }

		const int32_t us[4][13] =
		{
			{ 0, 1, 0, 1, 0,-1, 1,-1,-1, 2, 0, 2, 1},
			{ 0,-1, 0,-1, 0, 1,-1, 1, 1,-2, 0,-2,-1},
			{ 0, 0, 1, 1,-1, 0,-1, 1,-1, 0, 2, 1, 2},
			{ 0, 0,-1,-1, 1, 0, 1,-1, 1, 0,-2,-1,-2},
		};

		const int32_t vs[4][13] =
		{
			{ 0, 0, 1, 1,-1, 0,-1, 1,-1, 0, 2, 1, 2},
			{ 0, 0, 1, 1,-1, 0,-1, 1,-1, 0, 2, 1, 2},
			{ 0,-1, 0,-1, 0, 1,-1, 1, 1,-2, 0,-2,-1},
			{ 0,-1, 0,-1, 0, 1,-1, 1, 1,-2, 0,-2,-1},
		};

		auto sq_dist = std::numeric_limits<float>::max();
		auto cell_id = 0;

		// fast enough
		if (0.5f >= jitter)
		{
			for (auto loop = 0; loop < 4; ++loop)
			{
				const auto u       = us[quadrant][loop];
				const auto v       = vs[quadrant][loop];
				const auto hash    = my_hash(int_x + u, int_y + v);
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

		// cell id and possible minimum squared distance
		// 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12
		// 0, 0, 0, 0, 1, 1, 1, 1, 2, 4, 4, 4, 4
		const std::array<int32_t, 5> slices = { 0, 4, 8, 9, 13 };
		const std::array<int32_t, 4> ranges = { 0, 1, 2, 4 };

		const auto mine = my_hash(int_x, int_y);
		if (cache2d_cell_id == mine && cache2d_quadrant == quadrant)
		{
			// cache hit
			auto counter = 0;
			for (auto dist = 0; dist < 4; ++dist)
			{
				if (ranges[dist] < sq_dist * 4)
				{
					for (auto loop = slices[dist]; loop < slices[dist + 1]; ++loop)
					{
						if (counter < cache2d_counter)
						{
							// reuse cache
							const auto sample = cache2d_samples[counter];
							const auto tmp    = glm::dot(source - sample, source - sample);
							if (sq_dist > tmp)
							{
								sq_dist = tmp;
								cell_id = cache2d_cell_ids[counter];
							}
						}
						else
						{
							// fall back
							const auto u       = us[quadrant][loop];
							const auto v       = vs[quadrant][loop];
							const auto hash    = my_hash(int_x + u, int_y + v);
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
						++counter;
					}
				}
				else
				{
					break;
				}
			}
		}
		else
		{
			// cache miss
			cache2d_cell_id  = mine;
			cache2d_quadrant = quadrant;
			cache2d_counter  = 0;

			for (auto dist = 0; dist < 4; ++dist)
			{
				if (ranges[dist] < sq_dist * 4)
				{
					for (auto loop = slices[dist]; loop < slices[dist + 1]; ++loop)
					{
						const auto u       = us[quadrant][loop];
						const auto v       = vs[quadrant][loop];
						const auto hash    = my_hash(int_x + u, int_y + v);
						const auto offsetX = OffsetX(hash, jitter);
						const auto offsetY = OffsetY(hash, jitter);
						const auto sample  = origin + glm::vec2(u + offsetX, v + offsetY);
						const auto tmp     = glm::dot(source - sample, source - sample);

						// update cache
						cache2d_samples [cache2d_counter] = sample;
						cache2d_cell_ids[cache2d_counter] = hash;
						++cache2d_counter;

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
		}

		return std::make_pair(cell_id, sq_dist);
	}

	// naive honeycomb implementation
	static std::pair<int32_t, float> Evaluate2DHex(const glm::vec2& source, const int32_t(*my_hash)(const int32_t, const int32_t), const float jitter = 1.0f)
	{
		assert(0.0f <= jitter && jitter <= 1.0f);

		const auto one_3  = 1.0f / std::sqrt(3.0f);
		const auto local  = glm::vec2(glm::dot(glm::vec2(1.0f, -one_3), source), glm::dot(glm::vec2(0.0f, 2.0f * one_3), source));
		const auto flt_x  = std::floor(local.x);
		const auto flt_y  = std::floor(local.y);
		const auto int_x  = int32_t(flt_x);
		const auto int_y  = int32_t(flt_y);
		const auto origin = glm::vec2(flt_x, flt_y);

		auto sq_dist = std::numeric_limits<float>::max();
		auto cell_id = 0;

		// 1 (self) + 8 (neighbours)
		const auto size = 9;
		const std::array<int32_t, size> us = { 0, 0,-1, 1, 0,-1, 1,-1, 1 };
		const std::array<int32_t, size> vs = { 0,-1, 0, 0, 1,-1,-1, 1, 1 };

		// "A Low-Distortion Map Between Triangle and Square" by Eric Heitz
		// maps a unit - square point (x, y) to a unit - triangle point
		auto lower_triangle = [&](float& x, float& y)
		{
			if (y > x)
			{
				x *= 0.5f;
				y -= x;
			}
			else
			{
				y *= 0.5f;
				x -= y;
			}
		};

		auto upper_triangle = [&](float& x, float& y)
		{
			lower_triangle(x, y);
			x = 1 - x;
			y = 1 - x;
		};

		for (auto loop = 0; loop < size; ++loop)
		{
			const auto u    = us[loop];
			const auto v    = vs[loop];
			const auto hash = my_hash(int_x + u, int_y + v);

			// lower triangle
			{
				auto randomX = OffsetX(hash);
				auto randomY = OffsetY(hash);
				lower_triangle(randomX, randomY);
				const auto offsetX = randomX * jitter + 1.0f / 3.0f;
				const auto offsetY = randomY * jitter + 1.0f / 3.0f;
				const auto sample  = origin + glm::vec2(u + offsetX, v + offsetY);
				const auto global  = glm::vec2(glm::dot(glm::vec2(1.0f, 0.5f), sample), glm::dot(glm::vec2(0.0f, std::sqrt(3) * 0.5f), sample));
				const auto tmp     = glm::dot(source - global, source - global);
				if (sq_dist > tmp)
				{
					sq_dist = tmp;
					cell_id = hash;
				}
			}

			// upper triangle
			{
				auto randomX = OffsetX(hash);
				auto randomY = OffsetY(hash);
				upper_triangle(randomX, randomY);
				const auto offsetZ = randomX * jitter + 2.0f / 3.0f;
				const auto offsetT = randomY * jitter + 2.0f / 3.0f;
				const auto sample  = origin + glm::vec2(u + offsetZ, v + offsetT);
				const auto global  = glm::vec2(glm::dot(glm::vec2(), sample), glm::dot(glm::vec2(), sample));
				const auto tmp     = glm::dot(source - global, source - global);
				if (sq_dist > tmp)
				{
					sq_dist = tmp;
					cell_id = hash;
				}
			}
		}

		return std::make_pair(cell_id, sq_dist);
	}

	// TODO: optimized honeycomb implementation

	// naive implementation
	static std::pair<int, float> Evaluate3DRef(const glm::vec3& source, const int32_t (*my_hash)(const int32_t, const int32_t, const int32_t), const float jitter = 1.0f)
	{
		assert(0.0f <= jitter && jitter <= 1.0f);

		const auto flt_x  = std::floor(source.x);
		const auto flt_y  = std::floor(source.y);
		const auto flt_z  = std::floor(source.z);
		const auto int_x  = int32_t(flt_x);
		const auto int_y  = int32_t(flt_y);
		const auto int_z  = int32_t(flt_z);
		const auto origin = glm::vec3(flt_x, flt_y, flt_z);

		const auto size = 117;

		const std::array<int32_t, size> us = { 0, 0, 0,-1, 1, 0, 0, 0,-1, 1, 0,-1, 1,-1, 1, 0,-1, 1, 0,-1, 1,-1, 1,-1, 1,-1, 1, 0, 0,-2, 2, 0, 0, 0,-1, 1, 0, 0,-2, 2, 0,-1, 1,-2, 2,-2, 2,-1, 1, 0,-2, 2, 0, 0,-1, 1, 0,-1, 1,-1, 1,-1, 1,-2, 2,-2, 2,-1, 1,-1, 1,-2, 2,-2, 2,-1, 1,-1, 1,-1, 1, 0,-2, 2, 0,-2, 2,-2, 2, 0,-2, 2, 0,-1, 1,-2, 2,-2, 2,-1, 1,-2, 2,-2, 2,-2, 2,-2, 2,-1, 1,-2, 2,-2, 2,-1, 1};
		const std::array<int32_t, size> vs = { 0, 0,-1, 0, 0, 1, 0,-1, 0, 0, 1,-1,-1, 1, 1,-1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1, 1, 0,-2, 0, 0, 2, 0,-1, 0, 0, 1,-2, 0, 0, 2,-2,-2,-1,-1, 1, 1, 2, 2,-2, 0, 0, 2,-1, 0, 0, 1,-1,-1, 1, 1,-2,-2,-1,-1, 1, 1, 2, 2,-2,-2,-1,-1, 1, 1, 2, 2,-1,-1, 1, 1,-2, 0, 0, 2,-2,-2, 2, 2,-2, 0, 0, 2,-2,-2,-1,-1, 1, 1, 2, 2,-2,-2, 2, 2,-2,-2, 2, 2,-2,-2,-1,-1, 1, 1, 2, 2};
		const std::array<int32_t, size> ws = { 0,-1, 0, 0, 0, 0, 1,-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 1, 1, 1, 1,-2, 0, 0, 0, 0, 2,-2,-2,-2,-2,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,-2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,-2,-2,-2,-2, 0, 0, 0, 0, 2, 2, 2, 2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,-1,-1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2};

		auto sq_dist = std::numeric_limits<float>::max();
		auto cell_id = 0;
		for (auto loop = 0; loop < size; ++loop)
		{
			const auto u       = us[loop];
			const auto v       = vs[loop];
			const auto w       = ws[loop];
			const auto hash    = my_hash(int_x + u, int_y + v, int_z + w);
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
	static std::pair<int, float> Evaluate3DOpt(const glm::vec3& source, const int32_t (*my_hash)(const int32_t, const int32_t, const int32_t), const float jitter = 1.0f)
	{
		assert(0.0f <= jitter && jitter <= 1.0f);

		const auto flt_x  = std::floor(source.x);
		const auto flt_y  = std::floor(source.y);
		const auto flt_z  = std::floor(source.z);
		const auto int_x  = int32_t(flt_x);
		const auto int_y  = int32_t(flt_y);
		const auto int_z  = int32_t(flt_z);
		const auto origin = glm::vec3(flt_x, flt_y, flt_z);

		// early termination
		if (0.0f == jitter)
		{
			const auto hash   = my_hash(int_x, int_y, int_z);
			const auto sample = origin + 0.5f;
			return std::make_pair(hash, glm::dot(source - sample, source - sample));
		}

		// octant selection
		auto octant = 0;
		if (0.5f > source.x - origin.x) { octant += 1; }
		if (0.5f > source.y - origin.y) { octant += 2; }
		if (0.5f > source.z - origin.z) { octant += 4; }

		const int us[8][39] =
		{
			{ 0, 1, 0, 0, 1, 1, 0, 1, 0, 0,-1, 1, 0, 1,-1, 0,-1, 1, 1,-1, 0,-1,-1, 1,-1,-1,-1, 2, 0, 0, 2, 1, 2, 0, 1, 0, 2, 1, 1},
			{ 0,-1, 0, 0,-1,-1, 0,-1, 0, 0, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1, 0, 1, 1,-1, 1, 1, 1,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2,-1,-1},
			{ 0, 0, 1, 0, 1, 0, 1, 1, 0,-1, 0, 0, 1,-1, 1,-1, 0, 1,-1, 1,-1, 0,-1,-1, 1,-1,-1, 0, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2, 1},
			{ 0, 0,-1, 0,-1, 0,-1,-1, 0, 1, 0, 0,-1, 1,-1, 1, 0,-1, 1,-1, 1, 0, 1, 1,-1, 1, 1, 0,-2, 0,-1,-2, 0,-2, 0,-1,-1,-2,-1},
			{ 0, 0, 1, 0, 1, 0, 1, 1, 0,-1, 0, 0,-1, 1,-1, 1, 0, 1,-1, 1,-1, 0,-1,-1, 1,-1,-1, 0, 2, 0, 1, 0, 2, 0, 2, 1, 1, 2, 1},
			{ 0, 0,-1, 0,-1, 0,-1,-1, 0, 1, 0, 0, 1,-1, 1,-1, 0,-1, 1,-1, 1, 0, 1, 1,-1, 1, 1, 0,-2, 0,-1, 0,-2, 0,-2,-1,-1,-2,-1},
			{ 0, 0, 0, 1, 0, 1, 1, 1,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1,-1, 0,-1,-1, 1,-1, 0, 0, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2},
			{ 0, 0, 0,-1, 0,-1,-1,-1, 1, 0, 0, 1, 0, 1,-1, 0,-1, 1,-1,-1, 1, 1, 0, 1, 1,-1, 1, 0, 0,-2, 0,-1, 0,-2,-1,-2,-1,-1,-2},
		};

		const int vs[8][39] =
		{
			{ 0, 0, 1, 0, 1, 0, 1, 1, 0,-1, 0, 0, 1,-1, 1,-1, 0, 1,-1, 1,-1, 0,-1,-1, 1,-1,-1, 0, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2, 1},
			{ 0, 0, 1, 0, 1, 0, 1, 1, 0,-1, 0, 0, 1,-1, 1,-1, 0, 1,-1, 1,-1, 0,-1,-1, 1,-1,-1, 0, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2, 1},
			{ 0,-1, 0, 0,-1,-1, 0,-1, 0, 0, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1, 0, 1, 1,-1, 1, 1, 1,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2,-1,-1},
			{ 0,-1, 0, 0,-1,-1, 0,-1, 0, 0, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1, 0, 1, 1,-1, 1, 1, 1,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2,-1,-1},
			{ 0, 0, 0, 1, 0, 1, 1, 1,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1,-1, 0,-1,-1, 1,-1, 0, 0, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2},
			{ 0, 0, 0, 1, 0, 1, 1, 1,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1,-1, 0,-1,-1, 1,-1, 0, 0, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2},
			{ 0, 0,-1, 0,-1, 0,-1,-1, 0, 1, 0, 0, 1,-1, 1,-1, 0,-1, 1,-1, 1, 0, 1, 1,-1, 1, 1, 0,-2, 0,-1, 0,-2, 0,-2,-1,-1,-2,-1},
			{ 0, 0,-1, 0,-1, 0,-1,-1, 0, 1, 0, 0, 1,-1, 1,-1, 0,-1, 1,-1, 1, 0, 1, 1,-1, 1, 1, 0,-2, 0,-1, 0,-2, 0,-2,-1,-1,-2,-1},
		};

		const int ws[8][39] =
		{
			{ 0, 0, 0, 1, 0, 1, 1, 1,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1,-1, 0,-1,-1, 1,-1, 0, 0, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2},
			{ 0, 0, 0, 1, 0, 1, 1, 1,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1,-1, 0,-1,-1, 1,-1, 0, 0, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2},
			{ 0, 0, 0, 1, 0, 1, 1, 1,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1,-1, 0,-1,-1, 1,-1, 0, 0, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2},
			{ 0, 0, 0, 1, 0, 1, 1, 1,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1,-1, 0,-1,-1, 1,-1, 0, 0, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2},
			{ 0,-1, 0, 0,-1,-1, 0,-1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1, 0, 1, 1,-1, 1, 1, 1,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1},
			{ 0,-1, 0, 0,-1,-1, 0,-1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1, 0, 1, 1,-1, 1, 1, 1,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1},
			{ 0,-1, 0, 0,-1,-1, 0,-1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1, 0, 1, 1,-1, 1, 1, 1,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1},
			{ 0,-1, 0, 0,-1,-1, 0,-1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1, 0, 1, 1,-1, 1, 1, 1,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1},
		};

		auto sq_dist = std::numeric_limits<float>::max();
		auto cell_id = 0;

		if (0.5f >= jitter)
		{
			// 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19
			// 0, 1, 1, 1, 2, 2, 2, 3, 9, 9, 9,10,10,10,10,10,10,11,11,11
			const std::array<int32_t, 7> ranges = { 0, 1, 2, 3, 9, 10, 11 };
			const std::array<int32_t, 8> slices = { 0, 1, 4, 7, 8, 11, 17, 20 };

			for (auto dist = 0; dist < 7; ++dist)
			{
				if (ranges[dist] < sq_dist * 16)
				{
					for (auto loop = slices[dist]; loop < slices[dist + 1]; ++loop)
					{
						const auto u       = us[octant][loop];
						const auto v       = vs[octant][loop];
						const auto w       = ws[octant][loop];
						const auto hash    = my_hash(int_x + u, int_y + v, int_z + w);
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

		// cell id and possible minimum squared distance
		// 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38
		// 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
		const std::array<int32_t, 5> ranges = { 0,  1,  2,  3,  4 };
		const std::array<int32_t, 6> slices = { 0,  8, 20, 26, 27, 39 };

		for (auto dist = 0; dist < 5; ++dist)
		{
			if (ranges[dist] < sq_dist * 4)
			{
				for (auto loop = slices[dist]; loop < slices[dist + 1]; ++loop)
				{
					const auto u       = us[octant][loop];
					const auto v       = vs[octant][loop];
					const auto w       = ws[octant][loop];
					const auto hash    = my_hash(int_x + u, int_y + v, int_z + w);
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

	// naive implementation
	static std::pair<int, float> Evaluate4DRef(const glm::vec4& source, const int32_t (*my_hash)(const int32_t, const int32_t, const int32_t, const int32_t), const float jitter = 1.0f)
	{
		assert(0.0f <= jitter && jitter <= 1.0f);

		const auto flt_x  = std::floor(source.x);
		const auto flt_y  = std::floor(source.y);
		const auto flt_z  = std::floor(source.z);
		const auto flt_t  = std::floor(source.w);
		const auto int_x  = int32_t(flt_x);
		const auto int_y  = int32_t(flt_y);
		const auto int_z  = int32_t(flt_z);
		const auto int_t  = int32_t(flt_t);
		const auto origin = glm::vec4(flt_x, flt_y, flt_z, flt_t);

		auto sq_dist = std::numeric_limits<float>::max();
		auto cell_id = 0;
		for (auto t = -2; t <= 2; ++t)
		{
			for (auto w = -2; w <= 2; ++w)
			{
				for (auto v = -2; v <= 2; ++v)
				{
					for (auto u = -2; u <= 2; ++u)
					{
						if (2 == std::abs(u) && 2 == std::abs(v) && 2 == std::abs(w) && 2 == std::abs(t))
						{
							continue;
						}
						const auto hash    = my_hash(int_x + u, int_y + v, int_z + w, int_t + t);
						const auto offsetX = OffsetX(hash, jitter);
						const auto offsetY = OffsetY(hash, jitter);
						const auto offsetZ = OffsetZ(hash, jitter);
						const auto offsetT = OffsetT(hash, jitter);
						const auto sample  = origin + glm::vec4(u + offsetX, v + offsetY, w + offsetZ, t + offsetT);
						const auto tmp     = glm::dot(source - sample, source - sample);
						if (sq_dist > tmp)
						{
							sq_dist = tmp;
							cell_id = hash;
						}
					}
				}
			}
		}

		return std::make_pair(cell_id, sq_dist);
	}

	// hextant optimization
	static std::pair<int, float> Evaluate4DOpt(const glm::vec4& source, const int32_t(*my_hash)(const int32_t, const int32_t, const int32_t, const int32_t), const float jitter = 1.0f)
	{
		assert(0.0f <= jitter && jitter <= 1.0f);

		const auto flt_x  = std::floor(source.x);
		const auto flt_y  = std::floor(source.y);
		const auto flt_z  = std::floor(source.z);
		const auto flt_t  = std::floor(source.w);
		const auto int_x  = int32_t(flt_x);
		const auto int_y  = int32_t(flt_y);
		const auto int_z  = int32_t(flt_z);
		const auto int_t  = int32_t(flt_t);
		const auto origin = glm::vec4(flt_x, flt_y, flt_z, flt_t);

		auto sq_dist = std::numeric_limits<float>::max();
		auto cell_id = 0;

		// early termination
		if (0.0f == jitter)
		{
			const auto hash   = my_hash(int_x, int_y, int_z, int_t);
			const auto sample = origin + 0.5f;
			return std::make_pair(hash, glm::dot(source - sample, source - sample));
		}

		// hextant selection
		auto hextant = 0;
		if (0.5f > source.x - origin.x) { hextant += 1; }
		if (0.5f > source.y - origin.y) { hextant += 2; }
		if (0.5f > source.z - origin.z) { hextant += 4; }
		if (0.5f > source.t - origin.t) { hextant += 8; }

		const int us[16][195] =
		{
			{ 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0,-1, 0, 0,-1, 1, 0, 1,-1, 0,-1, 0, 0,-1, 1, 0, 1,-1, 0,-1, 1, 1,-1, 1, 0, 1,-1, 0,-1, 1, 1,-1, 1, 1,-1, 0, 0,-1, 0, 0,-1, 1, 0, 1,-1, 0,-1, 1, 0, 1,-1, 0,-1, 1, 1,-1, 1, 1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1, 2, 0, 0, 2, 0, 0, 2, 1, 2, 0, 1, 0, 2, 0, 0, 2, 1, 2, 0, 1, 0, 2, 1, 1, 2, 1, 2, 0, 1, 0, 2, 1, 1, 2, 1, 1, 2, 0, 0, 2, 0, 0, 2, 1, 2, 0, 1, 0, 2, 1, 2, 0, 1, 0, 2, 1, 1, 2, 1, 1},
			{ 0, 0,-1, 0, 0, 0,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1,-1, 0,-1,-1,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1, 0, 0,-1,-1, 0,-1,-1, 0,-1,-1, 0, 0, 1, 0, 0, 1,-1, 0,-1, 1, 0, 1, 0, 0, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1,-1,-1, 1, 0, 0, 1, 0, 0, 1,-1, 0,-1, 1, 0, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1,-1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 1,-2, 0, 0,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2,-1,-1,-2,-1,-2, 0,-1, 0,-2,-1,-1,-2,-1,-1,-2, 0, 0,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2,-1,-2, 0,-1, 0,-2,-1,-1,-2,-1,-1},
			{ 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0,-1, 0, 0,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 0, 0, 1,-1, 1,-1, 0, 1,-1, 1, 0, 1,-1, 1,-1, 0, 1,-1, 1, 1,-1, 1, 0,-1, 0, 0,-1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 1,-1, 1, 1,-1, 1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1,-1, 0, 2, 0, 0, 2, 0, 1, 2, 0, 2, 0, 1, 0, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2, 1, 1, 2, 0, 2, 0, 1, 1, 2, 1, 1, 2, 1, 0, 2, 0, 0, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2, 0, 2, 0, 1, 1, 2, 1, 1, 2, 1},
			{ 0, 0, 0,-1, 0, 0, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1,-1,-1, 0,-1, 0,-1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 0,-1,-1, 0,-1,-1,-1, 0, 1, 0, 0, 1, 0, 0,-1, 1,-1, 1, 0, 0, 1, 0, 0,-1, 1,-1, 1, 0,-1, 1,-1, 0,-1, 1,-1, 1, 0,-1, 1,-1,-1, 1,-1, 0, 1, 0, 0, 1, 0, 0,-1, 1,-1, 1, 0, 0,-1, 1,-1, 1, 0,-1, 1,-1,-1, 1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 0,-2, 0, 0,-2, 0,-1,-2, 0,-2, 0,-1, 0,-2, 0,-1,-2, 0,-2, 0,-1,-1,-2,-1,-1,-2, 0,-2, 0,-1,-1,-2,-1,-1,-2,-1, 0,-2, 0, 0,-2, 0,-1,-2, 0,-2, 0,-1,-1,-2, 0,-2, 0,-1,-1,-2,-1,-1,-2,-1},
			{ 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0,-1, 0, 0,-1, 0, 0,-1, 1,-1, 1, 0, 0,-1, 0, 0,-1, 1,-1, 1, 0, 1,-1, 1, 0,-1, 1,-1, 1, 0, 1,-1, 1, 1,-1, 1, 0,-1, 0, 0,-1, 0, 0,-1, 1,-1, 1, 0, 0,-1, 1,-1, 1, 0, 1,-1, 1, 1,-1, 1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1,-1, 0, 2, 0, 0, 2, 0, 1, 0, 2, 0, 2, 1, 0, 2, 0, 1, 0, 2, 0, 2, 1, 1, 2, 1, 1, 0, 2, 0, 2, 1, 1, 2, 1, 1, 2, 1, 0, 2, 0, 0, 2, 0, 1, 0, 2, 0, 2, 1, 1, 0, 2, 0, 2, 1, 1, 2, 1, 1, 2, 1},
			{ 0, 0, 0,-1, 0, 0, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1,-1,-1, 0,-1, 0,-1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 0,-1,-1, 0,-1,-1,-1, 0, 1, 0, 0, 1, 0, 0, 1,-1, 1,-1, 0, 0, 1, 0, 0, 1,-1, 1,-1, 0,-1, 1,-1, 0, 1,-1, 1,-1, 0,-1, 1,-1,-1, 1,-1, 0, 1, 0, 0, 1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0,-1, 1,-1,-1, 1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 0,-2, 0, 0,-2, 0,-1, 0,-2, 0,-2,-1, 0,-2, 0,-1, 0,-2, 0,-2,-1,-1,-2,-1,-1, 0,-2, 0,-2,-1,-1,-2,-1,-1,-2,-1, 0,-2, 0, 0,-2, 0,-1, 0,-2, 0,-2,-1,-1, 0,-2, 0,-2,-1,-1,-2,-1,-1,-2,-1},
			{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1, 1, 1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1, 1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1, 0, 0, 2, 0, 0, 2, 0, 1, 0, 2, 1, 2, 0, 0, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 1, 1, 2, 0, 0, 2, 0, 0, 2, 0, 1, 0, 2, 1, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 1, 1, 2},
			{ 0, 0, 0, 0,-1, 0, 0, 0,-1, 0,-1,-1, 0, 0,-1, 0,-1,-1,-1, 0,-1,-1, 0,-1,-1, 0, 0, 0,-1, 0, 0,-1, 0,-1,-1, 0,-1,-1,-1,-1, 1, 0, 0, 1, 0, 0, 1, 0, 1,-1, 0,-1, 1, 0, 0, 1, 0, 1,-1, 0,-1, 1,-1,-1, 1, 0, 1,-1, 0,-1, 1,-1,-1, 1,-1,-1, 1, 0, 0, 1, 0, 0, 1, 0, 1,-1, 0,-1, 1, 0, 1,-1, 0,-1, 1,-1,-1, 1,-1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 0, 0,-2, 0, 0,-2, 0,-1, 0,-2,-1,-2, 0, 0,-2, 0,-1, 0,-2,-1,-2,-1,-1,-2, 0,-1, 0,-2,-1,-2,-1,-1,-2,-1,-1,-2, 0, 0,-2, 0, 0,-2, 0,-1, 0,-2,-1,-2, 0,-1, 0,-2,-1,-2,-1,-1,-2,-1,-1,-2},
			{ 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0,-1, 0, 0,-1, 1, 0, 1,-1, 0,-1, 0, 0,-1, 1, 0, 1,-1, 0,-1, 1, 1,-1, 1, 0, 1,-1, 0,-1, 1, 1,-1, 1, 1,-1, 0, 0,-1, 0, 0,-1, 1, 0, 1,-1, 0,-1, 1, 0, 1,-1, 0,-1, 1, 1,-1, 1, 1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1,-1,-1, 2, 0, 0, 2, 0, 0, 2, 1, 2, 0, 1, 0, 2, 0, 0, 2, 1, 2, 0, 1, 0, 2, 1, 1, 2, 1, 2, 0, 1, 0, 2, 1, 1, 2, 1, 1, 2, 0, 0, 2, 0, 0, 2, 1, 2, 0, 1, 0, 2, 1, 2, 0, 1, 0, 2, 1, 1, 2, 1, 1},
			{ 0, 0,-1, 0, 0, 0,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1,-1, 0,-1,-1,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1, 0, 0,-1,-1, 0,-1,-1, 0,-1,-1, 0, 0, 1, 0, 0, 1,-1, 0,-1, 1, 0, 1, 0, 0, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1,-1,-1, 1, 0, 0, 1, 0, 0, 1,-1, 0,-1, 1, 0, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1,-1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 1,-2, 0, 0,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2,-1,-1,-2,-1,-2, 0,-1, 0,-2,-1,-1,-2,-1,-1,-2, 0, 0,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2,-1,-2, 0,-1, 0,-2,-1,-1,-2,-1,-1},
			{ 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0,-1, 0, 0,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 0, 0, 1,-1, 1,-1, 0, 1,-1, 1, 0, 1,-1, 1,-1, 0, 1,-1, 1, 1,-1, 1, 0,-1, 0, 0,-1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 1,-1, 1, 1,-1, 1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1,-1, 0, 2, 0, 0, 2, 0, 1, 2, 0, 2, 0, 1, 0, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2, 1, 1, 2, 0, 2, 0, 1, 1, 2, 1, 1, 2, 1, 0, 2, 0, 0, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2, 0, 2, 0, 1, 1, 2, 1, 1, 2, 1},
			{ 0, 0, 0,-1, 0, 0, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1,-1,-1, 0,-1, 0,-1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 0,-1,-1, 0,-1,-1,-1, 0, 1, 0, 0, 1, 0, 0,-1, 1,-1, 1, 0, 0, 1, 0, 0,-1, 1,-1, 1, 0,-1, 1,-1, 0,-1, 1,-1, 1, 0,-1, 1,-1,-1, 1,-1, 0, 1, 0, 0, 1, 0, 0,-1, 1,-1, 1, 0, 0,-1, 1,-1, 1, 0,-1, 1,-1,-1, 1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 0,-2, 0, 0,-2, 0,-1,-2, 0,-2, 0,-1, 0,-2, 0,-1,-2, 0,-2, 0,-1,-1,-2,-1,-1,-2, 0,-2, 0,-1,-1,-2,-1,-1,-2,-1, 0,-2, 0, 0,-2, 0,-1,-2, 0,-2, 0,-1,-1,-2, 0,-2, 0,-1,-1,-2,-1,-1,-2,-1},
			{ 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0,-1, 0, 0,-1, 0, 0,-1, 1,-1, 1, 0, 0,-1, 0, 0,-1, 1,-1, 1, 0, 1,-1, 1, 0,-1, 1,-1, 1, 0, 1,-1, 1, 1,-1, 1, 0,-1, 0, 0,-1, 0, 0,-1, 1,-1, 1, 0, 0,-1, 1,-1, 1, 0, 1,-1, 1, 1,-1, 1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1,-1, 0, 2, 0, 0, 2, 0, 1, 0, 2, 0, 2, 1, 0, 2, 0, 1, 0, 2, 0, 2, 1, 1, 2, 1, 1, 0, 2, 0, 2, 1, 1, 2, 1, 1, 2, 1, 0, 2, 0, 0, 2, 0, 1, 0, 2, 0, 2, 1, 1, 0, 2, 0, 2, 1, 1, 2, 1, 1, 2, 1},
			{ 0, 0, 0,-1, 0, 0, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1,-1,-1, 0,-1, 0,-1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 0,-1,-1, 0,-1,-1,-1, 0, 1, 0, 0, 1, 0, 0, 1,-1, 1,-1, 0, 0, 1, 0, 0, 1,-1, 1,-1, 0,-1, 1,-1, 0, 1,-1, 1,-1, 0,-1, 1,-1,-1, 1,-1, 0, 1, 0, 0, 1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0,-1, 1,-1,-1, 1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 0,-2, 0, 0,-2, 0,-1, 0,-2, 0,-2,-1, 0,-2, 0,-1, 0,-2, 0,-2,-1,-1,-2,-1,-1, 0,-2, 0,-2,-1,-1,-2,-1,-1,-2,-1, 0,-2, 0, 0,-2, 0,-1, 0,-2, 0,-2,-1,-1, 0,-2, 0,-2,-1,-1,-2,-1,-1,-2,-1},
			{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1, 1, 1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1, 1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1, 0, 0, 2, 0, 0, 2, 0, 1, 0, 2, 1, 2, 0, 0, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 1, 1, 2, 0, 0, 2, 0, 0, 2, 0, 1, 0, 2, 1, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 1, 1, 2},
			{ 0, 0, 0, 0,-1, 0, 0, 0,-1, 0,-1,-1, 0, 0,-1, 0,-1,-1,-1, 0,-1,-1, 0,-1,-1, 0, 0, 0,-1, 0, 0,-1, 0,-1,-1, 0,-1,-1,-1,-1, 1, 0, 0, 1, 0, 0, 1, 0, 1,-1, 0,-1, 1, 0, 0, 1, 0, 1,-1, 0,-1, 1,-1,-1, 1, 0, 1,-1, 0,-1, 1,-1,-1, 1,-1,-1, 1, 0, 0, 1, 0, 0, 1, 0, 1,-1, 0,-1, 1, 0, 1,-1, 0,-1, 1,-1,-1, 1,-1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 0, 0,-2, 0, 0,-2, 0,-1, 0,-2,-1,-2, 0, 0,-2, 0,-1, 0,-2,-1,-2,-1,-1,-2, 0,-1, 0,-2,-1,-2,-1,-1,-2,-1,-1,-2, 0, 0,-2, 0, 0,-2, 0,-1, 0,-2,-1,-2, 0,-1, 0,-2,-1,-2,-1,-1,-2,-1,-1,-2},
		};

		const int vs[16][195] =
		{
			{ 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0,-1, 0, 0,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 0, 0, 1,-1, 1,-1, 0, 1,-1, 1, 0, 1,-1, 1,-1, 0, 1,-1, 1, 1,-1, 1, 0,-1, 0, 0,-1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 1,-1, 1, 1,-1, 1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1,-1, 0, 2, 0, 0, 2, 0, 1, 2, 0, 2, 0, 1, 0, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2, 1, 1, 2, 0, 2, 0, 1, 1, 2, 1, 1, 2, 1, 0, 2, 0, 0, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2, 0, 2, 0, 1, 1, 2, 1, 1, 2, 1},
			{ 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0,-1, 0, 0,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 0, 0, 1,-1, 1,-1, 0, 1,-1, 1, 0, 1,-1, 1,-1, 0, 1,-1, 1, 1,-1, 1, 0,-1, 0, 0,-1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 1,-1, 1, 1,-1, 1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1,-1, 0, 2, 0, 0, 2, 0, 1, 2, 0, 2, 0, 1, 0, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2, 1, 1, 2, 0, 2, 0, 1, 1, 2, 1, 1, 2, 1, 0, 2, 0, 0, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2, 0, 2, 0, 1, 1, 2, 1, 1, 2, 1},
			{ 0, 0,-1, 0, 0, 0,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1,-1, 0,-1,-1,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1, 0, 0,-1,-1, 0,-1,-1, 0,-1,-1, 0, 0, 1, 0, 0, 1,-1, 0,-1, 1, 0, 1, 0, 0, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1,-1,-1, 1, 0, 0, 1, 0, 0, 1,-1, 0,-1, 1, 0, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1,-1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 1,-2, 0, 0,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2,-1,-1,-2,-1,-2, 0,-1, 0,-2,-1,-1,-2,-1,-1,-2, 0, 0,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2,-1,-2, 0,-1, 0,-2,-1,-1,-2,-1,-1},
			{ 0, 0,-1, 0, 0, 0,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1,-1, 0,-1,-1,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1, 0, 0,-1,-1, 0,-1,-1, 0,-1,-1, 0, 0, 1, 0, 0, 1,-1, 0,-1, 1, 0, 1, 0, 0, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1,-1,-1, 1, 0, 0, 1, 0, 0, 1,-1, 0,-1, 1, 0, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1,-1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 1,-2, 0, 0,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2,-1,-1,-2,-1,-2, 0,-1, 0,-2,-1,-1,-2,-1,-1,-2, 0, 0,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2,-1,-2, 0,-1, 0,-2,-1,-1,-2,-1,-1},
			{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1, 1, 1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1, 1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1, 0, 0, 2, 0, 0, 2, 0, 1, 0, 2, 1, 2, 0, 0, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 1, 1, 2, 0, 0, 2, 0, 0, 2, 0, 1, 0, 2, 1, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 1, 1, 2},
			{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1, 1, 1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1, 1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1, 0, 0, 2, 0, 0, 2, 0, 1, 0, 2, 1, 2, 0, 0, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 1, 1, 2, 0, 0, 2, 0, 0, 2, 0, 1, 0, 2, 1, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 1, 1, 2},
			{ 0, 0, 0,-1, 0, 0, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1,-1,-1, 0,-1, 0,-1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 0,-1,-1, 0,-1,-1,-1, 0, 1, 0, 0, 1, 0, 0, 1,-1, 1,-1, 0, 0, 1, 0, 0, 1,-1, 1,-1, 0,-1, 1,-1, 0, 1,-1, 1,-1, 0,-1, 1,-1,-1, 1,-1, 0, 1, 0, 0, 1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0,-1, 1,-1,-1, 1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 0,-2, 0, 0,-2, 0,-1, 0,-2, 0,-2,-1, 0,-2, 0,-1, 0,-2, 0,-2,-1,-1,-2,-1,-1, 0,-2, 0,-2,-1,-1,-2,-1,-1,-2,-1, 0,-2, 0, 0,-2, 0,-1, 0,-2, 0,-2,-1,-1, 0,-2, 0,-2,-1,-1,-2,-1,-1,-2,-1},
			{ 0, 0, 0,-1, 0, 0, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1,-1,-1, 0,-1, 0,-1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 0,-1,-1, 0,-1,-1,-1, 0, 1, 0, 0, 1, 0, 0, 1,-1, 1,-1, 0, 0, 1, 0, 0, 1,-1, 1,-1, 0,-1, 1,-1, 0, 1,-1, 1,-1, 0,-1, 1,-1,-1, 1,-1, 0, 1, 0, 0, 1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0,-1, 1,-1,-1, 1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 0,-2, 0, 0,-2, 0,-1, 0,-2, 0,-2,-1, 0,-2, 0,-1, 0,-2, 0,-2,-1,-1,-2,-1,-1, 0,-2, 0,-2,-1,-1,-2,-1,-1,-2,-1, 0,-2, 0, 0,-2, 0,-1, 0,-2, 0,-2,-1,-1, 0,-2, 0,-2,-1,-1,-2,-1,-1,-2,-1},
			{ 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0,-1, 0, 0,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 0, 0, 1,-1, 1,-1, 0, 1,-1, 1, 0, 1,-1, 1,-1, 0, 1,-1, 1, 1,-1, 1, 0,-1, 0, 0,-1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 1,-1, 1, 1,-1, 1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1,-1, 0, 2, 0, 0, 2, 0, 1, 2, 0, 2, 0, 1, 0, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2, 1, 1, 2, 0, 2, 0, 1, 1, 2, 1, 1, 2, 1, 0, 2, 0, 0, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2, 0, 2, 0, 1, 1, 2, 1, 1, 2, 1},
			{ 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0,-1, 0, 0,-1, 0, 0, 1,-1, 1,-1, 0, 0,-1, 0, 0, 1,-1, 1,-1, 0, 1,-1, 1, 0, 1,-1, 1,-1, 0, 1,-1, 1, 1,-1, 1, 0,-1, 0, 0,-1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 1,-1, 1, 1,-1, 1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1,-1, 0, 2, 0, 0, 2, 0, 1, 2, 0, 2, 0, 1, 0, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2, 1, 1, 2, 0, 2, 0, 1, 1, 2, 1, 1, 2, 1, 0, 2, 0, 0, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2, 0, 2, 0, 1, 1, 2, 1, 1, 2, 1},
			{ 0, 0,-1, 0, 0, 0,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1,-1, 0,-1,-1,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1, 0, 0,-1,-1, 0,-1,-1, 0,-1,-1, 0, 0, 1, 0, 0, 1,-1, 0,-1, 1, 0, 1, 0, 0, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1,-1,-1, 1, 0, 0, 1, 0, 0, 1,-1, 0,-1, 1, 0, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1,-1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 1,-2, 0, 0,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2,-1,-1,-2,-1,-2, 0,-1, 0,-2,-1,-1,-2,-1,-1,-2, 0, 0,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2,-1,-2, 0,-1, 0,-2,-1,-1,-2,-1,-1},
			{ 0, 0,-1, 0, 0, 0,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1,-1, 0,-1,-1,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1, 0, 0,-1,-1, 0,-1,-1, 0,-1,-1, 0, 0, 1, 0, 0, 1,-1, 0,-1, 1, 0, 1, 0, 0, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1,-1,-1, 1, 0, 0, 1, 0, 0, 1,-1, 0,-1, 1, 0, 1,-1, 0,-1, 1, 0, 1,-1,-1, 1,-1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 1,-2, 0, 0,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2,-1,-1,-2,-1,-2, 0,-1, 0,-2,-1,-1,-2,-1,-1,-2, 0, 0,-2, 0, 0,-2,-1,-2, 0,-1, 0,-2,-1,-2, 0,-1, 0,-2,-1,-1,-2,-1,-1},
			{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1, 1, 1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1, 1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1, 0, 0, 2, 0, 0, 2, 0, 1, 0, 2, 1, 2, 0, 0, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 1, 1, 2, 0, 0, 2, 0, 0, 2, 0, 1, 0, 2, 1, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 1, 1, 2},
			{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1, 1, 1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 1, 0, 1,-1, 0,-1, 1, 0, 1,-1, 1, 1,-1, 1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1, 0, 0, 2, 0, 0, 2, 0, 1, 0, 2, 1, 2, 0, 0, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 1, 1, 2, 0, 0, 2, 0, 0, 2, 0, 1, 0, 2, 1, 2, 0, 1, 0, 2, 1, 2, 1, 1, 2, 1, 1, 2},
			{ 0, 0, 0,-1, 0, 0, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1,-1,-1, 0,-1, 0,-1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 0,-1,-1, 0,-1,-1,-1, 0, 1, 0, 0, 1, 0, 0, 1,-1, 1,-1, 0, 0, 1, 0, 0, 1,-1, 1,-1, 0,-1, 1,-1, 0, 1,-1, 1,-1, 0,-1, 1,-1,-1, 1,-1, 0, 1, 0, 0, 1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0,-1, 1,-1,-1, 1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 0,-2, 0, 0,-2, 0,-1, 0,-2, 0,-2,-1, 0,-2, 0,-1, 0,-2, 0,-2,-1,-1,-2,-1,-1, 0,-2, 0,-2,-1,-1,-2,-1,-1,-2,-1, 0,-2, 0, 0,-2, 0,-1, 0,-2, 0,-2,-1,-1, 0,-2, 0,-2,-1,-1,-2,-1,-1,-2,-1},
			{ 0, 0, 0,-1, 0, 0, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1, 0,-1,-1,-1, 0,-1, 0,-1,-1, 0, 0,-1, 0, 0,-1, 0,-1, 0,-1,-1, 0,-1,-1,-1, 0, 1, 0, 0, 1, 0, 0, 1,-1, 1,-1, 0, 0, 1, 0, 0, 1,-1, 1,-1, 0,-1, 1,-1, 0, 1,-1, 1,-1, 0,-1, 1,-1,-1, 1,-1, 0, 1, 0, 0, 1, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0,-1, 1,-1,-1, 1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 0,-2, 0, 0,-2, 0,-1, 0,-2, 0,-2,-1, 0,-2, 0,-1, 0,-2, 0,-2,-1,-1,-2,-1,-1, 0,-2, 0,-2,-1,-1,-2,-1,-1,-2,-1, 0,-2, 0, 0,-2, 0,-1, 0,-2, 0,-2,-1,-1, 0,-2, 0,-2,-1,-1,-2,-1,-1,-2,-1},
		};

		const int ws[16][195] =
		{
			{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1,-1, 0, 0,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1, 1, 1,-1, 0, 0,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1, 1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1, 0, 0, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2, 0, 0, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2},
			{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1,-1, 0, 0,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1, 1, 1,-1, 0, 0,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1, 1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1, 0, 0, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2, 0, 0, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2},
			{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1,-1, 0, 0,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1, 1, 1,-1, 0, 0,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1, 1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1, 0, 0, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2, 0, 0, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2},
			{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1,-1, 0, 0,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1, 1, 1,-1, 0, 0,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1, 1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1, 0, 0, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2, 0, 0, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2},
			{ 0, 0,-1, 0, 0, 0,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1,-1, 0,-1,-1,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1, 0, 0,-1,-1, 0,-1,-1, 0,-1,-1, 0, 0, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 1, 0, 0, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 1,-2, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-1,-1,-2, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-1,-1},
			{ 0, 0,-1, 0, 0, 0,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1,-1, 0,-1,-1,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1, 0, 0,-1,-1, 0,-1,-1, 0,-1,-1, 0, 0, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 1, 0, 0, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 1,-2, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-1,-1,-2, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-1,-1},
			{ 0, 0,-1, 0, 0, 0,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1,-1, 0,-1,-1,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1, 0, 0,-1,-1, 0,-1,-1, 0,-1,-1, 0, 0, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 1, 0, 0, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 1,-2, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-1,-1,-2, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-1,-1},
			{ 0, 0,-1, 0, 0, 0,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1,-1, 0,-1,-1,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1, 0, 0,-1,-1, 0,-1,-1, 0,-1,-1, 0, 0, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 1, 0, 0, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 1,-2, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-1,-1,-2, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-1,-1},
			{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1,-1, 0, 0,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1, 1, 1,-1, 0, 0,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1, 1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1, 0, 0, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2, 0, 0, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2},
			{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1,-1, 0, 0,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1, 1, 1,-1, 0, 0,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1, 1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1, 0, 0, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2, 0, 0, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2},
			{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1,-1, 0, 0,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1, 1, 1,-1, 0, 0,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1, 1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1, 0, 0, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2, 0, 0, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2},
			{ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1,-1, 0, 0,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1, 1, 1,-1, 0, 0,-1, 0, 0,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1, 1, 1,-1, 1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1, 0,-1,-1, 0,-1,-1, 1,-1,-1, 1,-1,-1,-1,-1,-1, 0, 0, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2, 0, 0, 2, 0, 0, 2, 0, 0, 1, 1, 2, 2, 0, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2},
			{ 0, 0,-1, 0, 0, 0,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1,-1, 0,-1,-1,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1, 0, 0,-1,-1, 0,-1,-1, 0,-1,-1, 0, 0, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 1, 0, 0, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 1,-2, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-1,-1,-2, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-1,-1},
			{ 0, 0,-1, 0, 0, 0,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1,-1, 0,-1,-1,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1, 0, 0,-1,-1, 0,-1,-1, 0,-1,-1, 0, 0, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 1, 0, 0, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 1,-2, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-1,-1,-2, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-1,-1},
			{ 0, 0,-1, 0, 0, 0,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1,-1, 0,-1,-1,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1, 0, 0,-1,-1, 0,-1,-1, 0,-1,-1, 0, 0, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 1, 0, 0, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 1,-2, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-1,-1,-2, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-1,-1},
			{ 0, 0,-1, 0, 0, 0,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1,-1, 0,-1,-1,-1, 0, 0,-1,-1, 0,-1, 0, 0,-1, 0, 0,-1,-1, 0,-1,-1, 0,-1,-1, 0, 0, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 1, 0, 0, 1, 0, 0, 1,-1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 1, 1,-1,-1, 1,-1,-1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 0, 1, 1, 0, 1, 1,-1, 1, 1,-1, 1, 1, 1, 1, 1, 1, 1,-2, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-1,-1,-2, 0, 0,-2, 0, 0,-2,-2,-1,-1, 0, 0,-2,-2,-1,-1, 0, 0,-2,-1,-1,-2,-1,-1},
		};

		const int ts[16][195] =
		{
			{ 0,-1, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 1, 1, 1,-2,-1, 1, 2,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0,-1, 1,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2},
			{ 0,-1, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 1, 1, 1,-2,-1, 1, 2,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0,-1, 1,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2},
			{ 0,-1, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 1, 1, 1,-2,-1, 1, 2,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0,-1, 1,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2},
			{ 0,-1, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 1, 1, 1,-2,-1, 1, 2,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0,-1, 1,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2},
			{ 0,-1, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 1, 1, 1,-2,-1, 1, 2,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0,-1, 1,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2},
			{ 0,-1, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 1, 1, 1,-2,-1, 1, 2,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0,-1, 1,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2},
			{ 0,-1, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 1, 1, 1,-2,-1, 1, 2,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0,-1, 1,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2},
			{ 0,-1, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 1, 1, 1,-2,-1, 1, 2,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0,-1, 1,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2},
			{ 0,-1, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 1, 1, 1,-2,-1, 1, 2,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0,-1, 1,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2},
			{ 0,-1, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 1, 1, 1,-2,-1, 1, 2,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0,-1, 1,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2},
			{ 0,-1, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 1, 1, 1,-2,-1, 1, 2,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0,-1, 1,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2},
			{ 0,-1, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 1, 1, 1,-2,-1, 1, 2,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0,-1, 1,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2},
			{ 0,-1, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 1, 1, 1,-2,-1, 1, 2,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0,-1, 1,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2},
			{ 0,-1, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 1, 1, 1,-2,-1, 1, 2,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0,-1, 1,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2},
			{ 0,-1, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 1, 1, 1,-2,-1, 1, 2,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0,-1, 1,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2},
			{ 0,-1, 0, 0, 0, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 1, 1, 1,-2,-1, 1, 2,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2, 2, 2, 2, 0,-1, 1,-2, 2, 0, 0, 0,-1,-1,-1, 0, 0, 0, 0, 0, 0, 1, 1, 1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 1, 1, 1,-1,-1,-1, 1, 1, 1,-2,-2,-2, 2, 2, 2,-2,-2,-2,-2,-2,-2, 2, 2, 2, 2, 2, 2,-2,-2,-2, 2, 2, 2},
		};

		if (0.5f >= jitter)
		{
			// cell id and possible minimum squared distance
			//  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84
			//  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10
			const std::array<int32_t, 6> ranges = { 0, 1,  2,  3,  9, 10 };
			const std::array<int32_t, 7> slices = { 0, 5, 20, 25, 40, 55, 85 };

			// TODO: prove why 6 is ok
			for (auto dist = 0; dist < 6; ++dist)
			{
				if (ranges[dist] < sq_dist * 16)
				{
					for (auto loop = slices[dist]; loop < slices[dist + 1]; ++loop)
					{
						const auto u       = us[hextant][loop];
						const auto v       = vs[hextant][loop];
						const auto w       = ws[hextant][loop];
						const auto t       = ts[hextant][loop];
						const auto hash    = my_hash(int_x + u, int_y + v, int_z + w, int_t + t);
						const auto offsetX = OffsetX(hash, jitter);
						const auto offsetY = OffsetY(hash, jitter);
						const auto offsetZ = OffsetZ(hash, jitter);
						const auto offsetT = OffsetT(hash, jitter);
						const auto sample  = origin + glm::vec4(u + offsetX, v + offsetY, w + offsetZ, t + offsetT);
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

		// cell id and possible minimum squared distance
		// 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194
		// 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4
		const std::array<int32_t, 5> ranges = { 0,   1,   2,   3,   4 };
		const std::array<int32_t, 6> slices = { 0,  40, 100, 130, 135, 195 };

		for (auto dist = 0; dist < 5; ++dist)
		{
			if (ranges[dist] < sq_dist * 4)
			{
				for (auto loop = slices[dist]; loop < slices[dist + 1]; ++loop)
				{
					const auto u       = us[hextant][loop];
					const auto v       = vs[hextant][loop];
					const auto w       = ws[hextant][loop];
					const auto t       = ts[hextant][loop];
					const auto hash    = my_hash(int_x + u, int_y + v, int_z + w, int_t + t);
					const auto offsetX = OffsetX(hash, jitter);
					const auto offsetY = OffsetY(hash, jitter);
					const auto offsetZ = OffsetZ(hash, jitter);
					const auto offsetT = OffsetT(hash, jitter);
					const auto sample  = origin + glm::vec4(u + offsetX, v + offsetY, w + offsetZ, t + offsetT);
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
