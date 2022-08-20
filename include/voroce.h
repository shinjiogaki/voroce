// APACHE LICENSE, VERSION 2.0
// copyright (c) shinji ogaki

#ifndef VOROCE
#define VOROCE

#include "../glm/glm/glm.hpp"

#include <algorithm>
#include <array>
#include <tuple>

namespace voroce
{
struct Voronoi
{
	// hash function from "Optimized Spatial Hashing for Collision Detection of Deformable Objects"
	// https://matthias-research.github.io/pages/publications/tetraederCollision.pdf
	static const auto PrimeU = 73856093;
	static const auto PrimeV = 13467053;
	static const auto PrimeW = 83492791;
	static const auto PrimeT = 23469181;
	static const auto LCG = 48271;
	static int32_t Hash2DLowQuality(const glm::ivec2& p)
	{
		return ((p.x * PrimeU) ^ (p.y * PrimeV)) * LCG;
	}
	static int32_t Hash3DLowQuality(const glm::ivec3& p)
	{
		return ((p.x * PrimeU) ^ (p.y * PrimeV) ^ (p.z * PrimeW)) * LCG;
	}
	static int32_t Hash4DLowQuality(const glm::ivec4& p)
	{
		return ((p.x * PrimeU) ^ (p.y * PrimeV) ^ (p.z * PrimeW) ^ (p.w * PrimeT)) * LCG;
	}

	// reference implementation
	static std::tuple<int32_t, float, glm::vec2> Evaluate2DRef(const glm::vec2& source, int32_t(*my_hash)(const glm::ivec2& p), const float jitter = 1.0f);
	static std::tuple<int32_t, float, glm::vec3> Evaluate3DRef(const glm::vec3& source, int32_t(*my_hash)(const glm::ivec3& p), const float jitter = 1.0f);
	static std::tuple<int32_t, float, glm::vec4> Evaluate4DRef(const glm::vec4& source, int32_t(*my_hash)(const glm::ivec4& p), const float jitter = 1.0f);

	// optimized (for bounced rays)
	static std::tuple<int32_t, float, glm::vec2> Evaluate2DOpt(const glm::vec2& source, int32_t(*my_hash)(const glm::ivec2& p), const float jitter = 1.0f);
	static std::tuple<int32_t, float, glm::vec3> Evaluate3DOpt(const glm::vec3& source, int32_t(*my_hash)(const glm::ivec3& p), const float jitter = 1.0f);
	static std::tuple<int32_t, float, glm::vec4> Evaluate4DOpt(const glm::vec4& source, int32_t(*my_hash)(const glm::ivec4& p), const float jitter = 1.0f);

	// additional optimization with cach (dedicated for primary rays)
	std::tuple<int32_t, float, glm::vec2> Evaluate2DCache(const glm::vec2& source, int32_t(*my_hash)(const glm::ivec2& p), const float jitter = 1.0f);
	std::tuple<int32_t, float, glm::vec3> Evaluate3DCache(const glm::vec3& source, int32_t(*my_hash)(const glm::ivec3& p), const float jitter = 1.0f);
	std::tuple<int32_t, float, glm::vec4> Evaluate4DCache(const glm::vec4& source, int32_t(*my_hash)(const glm::ivec4& p), const float jitter = 1.0f);

	// non rectangular grids
	static std::tuple<int32_t, float, glm::vec2> Evaluate2DTri(const glm::vec2& source, int32_t(*my_hash)(const glm::ivec2& p), const float jitter = 1.0f);
	static std::tuple<int32_t, float, glm::vec2> Evaluate2DHex(const glm::vec2& source, int32_t(*my_hash)(const glm::ivec2& p), const float jitter = 1.0f);

	// voronoi edge with normal vector
	// voronoi edge was originaly developed by Inigo Quilez https://iquilezles.org/articles/voronoilines/
	// voroce will use optimized neighbourhood query (TODO)
	// voroce provides normal vectors
	static std::tuple<int32_t, float, glm::vec2> Edge2DSharp(const glm::vec2& source, int32_t (*my_hash)(const glm::ivec2& p), const float jitter = 1.0f, const float width = 0.0f);
	static std::tuple<int32_t, float, glm::vec2> Edge2DRound(const glm::vec2& source, int32_t (*my_hash)(const glm::ivec2& p), const float jitter = 1.0f, const float width = 0.0f);

	// returns (cell_id, distance to the boundary, direction to the closest boundary)
	static auto Edge2DSharp(const float x, const float y, const float jitter = 1.0f, const float width = 0.0f) { return Voronoi::Edge2DSharp(glm::vec2(x, y), Voronoi::Hash2DLowQuality, jitter, width); }
	static auto Edge2DRound(const float x, const float y, const float jitter = 1.0f, const float width = 0.0f) { return Voronoi::Edge2DRound(glm::vec2(x, y), Voronoi::Hash2DLowQuality, jitter, width); }

private:
	// constants
	static const float pi;
	static const float sqrt3;
	static const float one_over_sqrt3;

	// minstd_rand (TODO: do something better, these are fast but a bit ugly...)
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

	// traversal order
	static const int32_t us2[4][13];
	static const int32_t vs2[4][13];
	static const int32_t us3[8][39];
	static const int32_t vs3[8][39];
	static const int32_t ws3[8][39];
	static const int32_t us4[16][195];
	static const int32_t vs4[16][195];
	static const int32_t ws4[16][195];
	static const int32_t ts4[16][195];

	// cache for 2d
	int32_t   cache2d_cell_id  = 0xFFFFFFFF;
	int32_t   cache2d_quadrant = 0;
	int32_t   cache2d_counter  = 0;
	int32_t   cache2d_cell_ids[13];
	glm::vec2 cache2d_samples [13];

	// cache for 3d
	int32_t   cache3d_cell_id = 0xFFFFFFFF;
	int32_t   cache3d_octant  = 0;
	int32_t   cache3d_counter = 0;
	int32_t   cache3d_cell_ids[39];
	glm::vec3 cache3d_samples [39];

	// cache for 4d
	int32_t   cache4d_cell_id = 0xFFFFFFFF;
	int32_t   cache4d_hextant = 0;
	int32_t   cache4d_counter = 0;
	int32_t   cache4d_cell_ids[195];
	glm::vec4 cache4d_samples [195];

	// SDF (line segment) by Inigo Quilez https://iquilezles.org/articles/distfunctions2d/
	static auto sd_segment(const glm::vec2 &p, const glm::vec2 &a, const glm::vec2 &b)
	{
		auto pa = p - a, ba = b - a;
		auto h = std::clamp(glm::dot(pa, ba) / glm::dot(ba, ba), 0.0f, 1.0f);
		return glm::length(pa - ba * h);
	}

	// SDF (Bezier curve) by Inigo Quilez https://iquilezles.org/articles/distfunctions2d/
	static auto sd_bezier(const glm::vec2 &pos, const glm::vec2 &A, const glm::vec2 &B, const glm::vec2 &C)
	{
		auto dot2  = [](const glm::vec2 &a) { return glm::dot(a, a); };
		auto cross = [](const glm::vec2 &a, const glm::vec2 &b) { return a.x * b.y - a.y * b.x; };

		auto a = B - A;
		auto b = A - 2.0f * B + C;
		auto c = a * 2.0f;
		auto d = A - pos;
		auto kk = 1.0f  / glm::dot(b, b);
		auto kx = kk * glm::dot(a, b);
		auto ky = kk * (2.0f * glm::dot(a, a) + glm::dot(d, b)) / 3.0f;
		auto kz = kk * glm::dot(d, a);
		auto p = ky - kx * kx;
		auto p3 = p * p * p;
		auto q = kx * (2.0f * kx * kx - 3.0f * ky) + kz;
		auto h = q * q + 4.0f * p3;
		if(0.0f <= h)
		{
			h = std::sqrt(h);
			auto x = (glm::vec2(h, -h) - q) / 2.0f;
			auto uv = glm::sign(x) * glm::pow(glm::abs(x), glm::vec2(1.0f / 3.0f));
			auto t  = glm::clamp(uv.x + uv.y - kx, 0.0f, 1.0f);
			auto q  = d + (c + b * t) * t;
			auto res = dot2(q);
			auto sgn = cross(c + 2.0f * b * t, q);
			return std::make_tuple(std::sqrt(res) * glm::sign(sgn), q);
		}
		auto z = std::sqrt(-p);
		auto v = std::acos(q / (p * z * 2.0f)) / 3.0f;
		auto m = std::cos(v);
		auto n = std::sin(v) * sqrt3;
		auto t = glm::clamp(glm::vec3(m + m, -n - m, n - m) * z - kx, 0.0f, 1.0f);
		auto qx = d + (c + b * t.x) * t.x;
		auto qy = d + (c + b * t.y) * t.y;
		auto dx = dot2(qx), sx = cross(c + 2.0f * b * t.x, qx);
		auto dy = dot2(qy), sy = cross(c + 2.0f * b * t.y, qy);
		if(dx < dy)
		{
			return std::make_tuple(std::sqrt(dx) * glm::sign(sx), qx);
		}
		return std::make_tuple(std::sqrt(dy) * glm::sign(sy), qy);
	}
};

#ifdef VOROCE_IMPLEMENTATION

const float Voronoi::pi             = 3.1415926535f;
const float Voronoi::sqrt3          = std::sqrt(3.0f);
const float Voronoi::one_over_sqrt3 = 1.0f / std::sqrt(3.0f);

// naive implementation
std::tuple<int32_t, float, glm::vec2> Voronoi::Evaluate2DRef(const glm::vec2& source, int32_t (*my_hash)(const glm::ivec2 &p), const float jitter)
{
	assert(0.0f <= jitter && jitter <= 1.0f);

	const auto origin    = glm::vec2(std::floor(source.x), std::floor(source.y));
	const auto quantized = glm::ivec2(int32_t(origin.x), int32_t(origin.y));

	auto sq_dist = std::numeric_limits<float>::max();
	auto cell_id = 0;
	glm::vec2 point;

	// 1 (self) + 20 (neighbours)
	const auto size = 21;
	const std::array<int32_t, size> us = { 0, 0,-1, 1, 0,-1, 1,-1, 1, 0,-2, 2, 0,-1, 1,-2, 2,-2, 2,-1, 1 };
	const std::array<int32_t, size> vs = { 0,-1, 0, 0, 1,-1,-1, 1, 1,-2, 0, 0, 2,-2,-2,-1,-1, 1, 1, 2, 2 };

	for (auto loop = 0; loop < size; ++loop)
	{
		const auto shift  = glm::ivec2(us[loop], vs[loop]);
		const auto hash   = my_hash(quantized + shift);
		const auto offset = glm::vec2(OffsetX(hash, jitter), OffsetY(hash, jitter));
		const auto sample = origin + offset + glm::vec2(shift);
		const auto tmp    = glm::dot(source - sample, source - sample);
		if (sq_dist > tmp)
		{
			sq_dist = tmp;
			cell_id = hash;
			point   = sample;
		}
	}

	return std::make_tuple(cell_id, sq_dist, point);
}

// quadrant optimization
std::tuple<int32_t, float, glm::vec2> Voronoi::Evaluate2DOpt(const glm::vec2& source, int32_t (*my_hash)(const glm::ivec2& p), const float jitter)
{
	assert(0.0f <= jitter && jitter <= 1.0f);

	const auto origin    = glm::vec2(std::floor(source.x), std::floor(source.y));
	const auto quantized = glm::ivec2(int32_t(origin.x), int32_t(origin.y));

	// early termination
	if (0.0f == jitter)
	{
		const auto hash   = my_hash(quantized);
		const auto sample = origin + 0.5f;
		return std::make_tuple(hash, glm::dot(source - sample, source - sample), sample);
	}

	// quadrant selection
	auto quadrant = 0;
	if (0.5f > source.x - origin.x) { quadrant += 1; }
	if (0.5f > source.y - origin.y) { quadrant += 2; }

	auto sq_dist = std::numeric_limits<float>::max();
	auto cell_id = 0;
	glm::vec2 point;

	if (0.5f >= jitter)
	{
		for (auto loop = 0; loop < 4; ++loop)
		{
			const auto shift  = glm::ivec2(us2[quadrant][loop], vs2[quadrant][loop]);
			const auto hash   = my_hash(quantized + shift);
			const auto offset = glm::vec2(OffsetX(hash, jitter), OffsetY(hash, jitter));
			const auto sample = origin + offset + glm::vec2(shift);
			const auto tmp    = glm::dot(source - sample, source - sample);
			if (sq_dist > tmp)
			{
				sq_dist = tmp;
				cell_id = hash;
				point   = sample;
			}
		}

		return std::make_tuple(cell_id, sq_dist, point);
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
				const auto shift  = glm::ivec2(us2[quadrant][loop], vs2[quadrant][loop]);
				const auto hash   = my_hash(quantized + shift);
				const auto offset = glm::vec2(OffsetX(hash, jitter), OffsetY(hash, jitter));
				const auto sample = origin + offset + glm::vec2(shift);
				const auto tmp    = glm::dot(source - sample, source - sample);
				if (sq_dist > tmp)
				{
					sq_dist = tmp;
					cell_id = hash;
					point   = sample;
				}
			}
		}
		else
		{
			break;
		}
	}

	return std::make_tuple(cell_id, sq_dist, point);
}

// quadrant optimization with cache
std::tuple<int32_t, float, glm::vec2> Voronoi::Evaluate2DCache(const glm::vec2& source, int32_t(*my_hash)(const glm::ivec2& p), const float jitter)
{
	assert(0.0f <= jitter && jitter <= 1.0f);

	const auto origin    = glm::vec2(std::floor(source.x), std::floor(source.y));
	const auto quantized = glm::ivec2(int32_t(origin.x), int32_t(origin.y));
	const auto self      = my_hash(quantized);

	// early termination
	if (0.0f == jitter)
	{
		const auto sample = origin + 0.5f;
		return std::make_tuple(self, glm::dot(source - sample, source - sample), sample);
	}

	// quadrant selection
	auto quadrant = 0;
	if (0.5f > source.x - origin.x) { quadrant += 1; }
	if (0.5f > source.y - origin.y) { quadrant += 2; }

	auto sq_dist = std::numeric_limits<float>::max();
	auto cell_id = 0;
	glm::vec2 point;

	// fast enough
	if (0.5f >= jitter)
	{
		if (cache2d_cell_id == self && cache2d_quadrant == quadrant)
		{
			// cache hit
			auto counter = 0;
			for (auto loop = 0; loop < 4; ++loop)
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
						point = sample;
					}
				}
				else
				{
					// fall back
					const auto shift  = glm::ivec2(us2[quadrant][loop], vs2[quadrant][loop]);
					const auto hash   = my_hash(quantized + shift);
					const auto offset = glm::vec2(OffsetX(hash, jitter), OffsetY(hash, jitter));
					const auto sample = origin + offset + glm::vec2(shift);
					const auto tmp    = glm::dot(source - sample, source - sample);
					if (sq_dist > tmp)
					{
						sq_dist = tmp;
						cell_id = hash;
						point   = sample;
					}
				}
				++counter;
			}
		}
		else
		{
			// cache miss
			cache2d_cell_id  = self;
			cache2d_quadrant = quadrant;
			cache2d_counter  = 0;

			for (auto loop = 0; loop < 4; ++loop)
			{
				const auto shift  = glm::ivec2(us2[quadrant][loop], vs2[quadrant][loop]);
				const auto hash   = my_hash(quantized + shift);
				const auto offset = glm::vec2(OffsetX(hash, jitter), OffsetY(hash, jitter));
				const auto sample = origin + offset + glm::vec2(shift);
				const auto tmp    = glm::dot(source - sample, source - sample);

				// update cache
				cache2d_samples [cache2d_counter] = sample;
				cache2d_cell_ids[cache2d_counter] = hash;
				++cache2d_counter;

				if (sq_dist > tmp)
				{
					sq_dist = tmp;
					cell_id = hash;
					point   = sample;
				}
			}
		}
		return std::make_tuple(cell_id, sq_dist, point);
	}

	// cell id and possible minimum squared distance
	// 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12
	// 0, 0, 0, 0, 1, 1, 1, 1, 2, 4, 4, 4, 4
	const std::array<int32_t, 5> slices = { 0, 4, 8, 9, 13 };
	const std::array<int32_t, 4> ranges = { 0, 1, 2, 4 };

	if (cache2d_cell_id == self && cache2d_quadrant == quadrant)
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
							point   = sample;
						}
					}
					else
					{
						// fall back
						const auto shift  = glm::ivec2(us2[quadrant][loop], vs2[quadrant][loop]);
						const auto hash   = my_hash(quantized + shift);
						const auto offset = glm::vec2(OffsetX(hash, jitter), OffsetY(hash, jitter));
						const auto sample = origin + offset + glm::vec2(shift);
						const auto tmp    = glm::dot(source - sample, source - sample);
						if (sq_dist > tmp)
						{
							sq_dist = tmp;
							cell_id = hash;
							point   = sample;
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
		cache2d_cell_id  = self;
		cache2d_quadrant = quadrant;
		cache2d_counter  = 0;

		for (auto dist = 0; dist < 4; ++dist)
		{
			if (ranges[dist] < sq_dist * 4)
			{
				for (auto loop = slices[dist]; loop < slices[dist + 1]; ++loop)
				{
					const auto shift  = glm::ivec2(us2[quadrant][loop], vs2[quadrant][loop]);
					const auto hash   = my_hash(quantized + shift);
					const auto offset = glm::vec2(OffsetX(hash, jitter), OffsetY(hash, jitter));
					const auto sample = origin + offset + glm::vec2(shift);
					const auto tmp    = glm::dot(source - sample, source - sample);

					// update cache
					cache2d_samples [cache2d_counter] = sample;
					cache2d_cell_ids[cache2d_counter] = hash;
					++cache2d_counter;

					if (sq_dist > tmp)
					{
						sq_dist = tmp;
						cell_id = hash;
						point   = sample;
					}
				}
			}
			else
			{
				break;
			}
		}
	}

	return std::make_tuple(cell_id, sq_dist, point);
}

// sharp Voronoi edges
std::tuple<int32_t, float, glm::vec2> Voronoi::Edge2DSharp(const glm::vec2& shading, int32_t (*my_hash)(const glm::ivec2& p), const float jitter, const float width)
{
	assert(0.0f <= jitter && jitter <= 1.0f);

	const auto [old_cell_id, sq_dis, closest] = Voronoi::Evaluate2DOpt(shading, my_hash, jitter);

	const auto origin    = glm::vec2(std::floor(closest.x), std::floor(closest.y));
	const auto quantized = glm::ivec2(int32_t(origin.x), int32_t(origin.y));

	// quadrant selection
	auto quadrant = 0;
	if (0.5f > closest.x - origin.x) { quadrant += 1; }
	if (0.5f > closest.y - origin.y) { quadrant += 2; }

	// list up all offsets
	auto dis = std::numeric_limits<float>::max();
	glm::vec2 point;
	const auto N = 13;
	for (auto loop = 1; loop < N; ++loop)
	{
		const auto shift  = glm::ivec2(us2[quadrant][loop], vs2[quadrant][loop]);
		const auto hash   = my_hash(quantized + shift);
		const auto offset = glm::vec2(OffsetX(hash, jitter), OffsetY(hash, jitter));
		const auto sample = origin + offset + glm::vec2(shift);
		const auto tmp    = glm::dot(shading - (sample + closest) * 0.5f, glm::normalize(closest - sample));
		if (dis > tmp)
		{
			dis   = tmp;
			point = sample;
		}
	}
	if(width > dis)
	{
		const auto dir = (point - closest) * 0.5f;
		return std::make_tuple(old_cell_id, dis, dir);
	}
	return std::make_tuple(old_cell_id, dis, glm::vec2(0));
}

// round Voronoi edges
std::tuple<int32_t, float, glm::vec2> Voronoi::Edge2DRound(const glm::vec2& shading, int32_t (*my_hash)(const glm::ivec2& p), const float jitter, const float width)
{
	assert(0.0f <= jitter && jitter <= 1.0f);

	const auto [old_cell_id, sq_dis, closest] = Voronoi::Evaluate2DOpt(shading, my_hash, jitter);

	const auto origin    = glm::vec2(std::floor(closest.x), std::floor(closest.y));
	const auto quantized = glm::ivec2(int32_t(origin.x), int32_t(origin.y));

	// quadrant selection
	auto quadrant = 0;
	if (0.5f > closest.x - origin.x) { quadrant += 1; }
	if (0.5f > closest.y - origin.y) { quadrant += 2; }

	// list up all offsets
	glm::vec2 point;
	const auto size = 13;
	const auto N = size - 1;
	std::array<glm::vec2, size> offsets;
	for (auto loop = 1; loop < size; ++loop)
	{
		const auto shift  = glm::ivec2(us2[quadrant][loop], vs2[quadrant][loop]);
		const auto hash   = my_hash(quantized + shift);
		const auto offset = glm::vec2(OffsetX(hash, jitter), OffsetY(hash, jitter));
		const auto sample = origin + offset + glm::vec2(shift);
		offsets[loop - 1] = sample - closest;
	}

	// compute intersections
	std::vector<glm::vec2> intersects; intersects.reserve(32);
	for (auto i = 0; i < N; ++i)
	{
		for (auto j = i + 1; j < N; ++j)
		{
			auto matrix = glm::mat2x2(offsets[i].x, offsets[i].y, offsets[j].x, offsets[j].y);
			if(0.0f == glm::determinant(matrix))
			{
				continue;
			}
			auto intersect = glm::vec2(glm::dot(offsets[i],offsets[i]), glm::dot(offsets[j], offsets[j])) * glm::inverse(matrix);
			auto flag      = true;
			for (auto k = 0; k < N; ++k) if (i != k && j != k)
			{
				if (0.0f < glm::dot(offsets[k], intersect - offsets[k]))
				{
					flag = false;
					break;
				}
			}
			if (flag)
			{
				intersects.push_back(intersect);
			}
		}
	}

	// sort intersections
	auto polar = [&](const glm::vec2 &p){ return (std::atan2(p.y, p.x) + pi) / (2 * pi); };
	std::sort(begin(intersects), end(intersects), [&](const glm::vec2& a, const glm::vec2& b) { return polar(a) > polar(b); });

	// bezier
	auto sgn = std::numeric_limits<float>::max();
	auto dis = std::numeric_limits<float>::max();
	auto dir = glm::vec2(0);
	auto I = int32_t(std::size(intersects));
	for (auto i = 0; i < I; ++i)
	{
		const auto p = (i - 1 + I) % I;
		const auto n = (i + 1    ) % I;
		const auto mid_p = (intersects[i] + intersects[p]) * 0.5f;
		const auto mid_n = (intersects[i] + intersects[n]) * 0.5f;
		const auto [d, direction] = sd_bezier(shading - closest, mid_p * 0.5f, intersects[i] * 0.5f, mid_n * 0.5f);
		if (sgn > std::abs(d))
		{
			dis = d;
			sgn = std::abs(d);
			dir = direction;
		}
	}
	if(width > dis)
	{
		return std::make_tuple(old_cell_id, dis, dir);
	}
	return std::make_tuple(old_cell_id, dis, glm::vec2(0));
}

// naive implementation
std::tuple<int32_t, float, glm::vec3> Voronoi::Evaluate3DRef(const glm::vec3& source, int32_t (*my_hash)(const glm::ivec3& p), const float jitter)
{
	assert(0.0f <= jitter && jitter <= 1.0f);

	const auto origin    = glm::vec3(std::floor(source.x), std::floor(source.y), std::floor(source.z));
	const auto quantized = glm::ivec3(int32_t(origin.x), int32_t(origin.y), int32_t(origin.z));

	const auto size = 117;
	const std::array<int32_t, size> us = { 0, 0, 0,-1, 1, 0, 0, 0,-1, 1, 0,-1, 1,-1, 1, 0,-1, 1, 0,-1, 1,-1, 1,-1, 1,-1, 1, 0, 0,-2, 2, 0, 0, 0,-1, 1, 0, 0,-2, 2, 0,-1, 1,-2, 2,-2, 2,-1, 1, 0,-2, 2, 0, 0,-1, 1, 0,-1, 1,-1, 1,-1, 1,-2, 2,-2, 2,-1, 1,-1, 1,-2, 2,-2, 2,-1, 1,-1, 1,-1, 1, 0,-2, 2, 0,-2, 2,-2, 2, 0,-2, 2, 0,-1, 1,-2, 2,-2, 2,-1, 1,-2, 2,-2, 2,-2, 2,-2, 2,-1, 1,-2, 2,-2, 2,-1, 1};
	const std::array<int32_t, size> vs = { 0, 0,-1, 0, 0, 1, 0,-1, 0, 0, 1,-1,-1, 1, 1,-1, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1, 1, 0,-2, 0, 0, 2, 0,-1, 0, 0, 1,-2, 0, 0, 2,-2,-2,-1,-1, 1, 1, 2, 2,-2, 0, 0, 2,-1, 0, 0, 1,-1,-1, 1, 1,-2,-2,-1,-1, 1, 1, 2, 2,-2,-2,-1,-1, 1, 1, 2, 2,-1,-1, 1, 1,-2, 0, 0, 2,-2,-2, 2, 2,-2, 0, 0, 2,-2,-2,-1,-1, 1, 1, 2, 2,-2,-2, 2, 2,-2,-2, 2, 2,-2,-2,-1,-1, 1, 1, 2, 2};
	const std::array<int32_t, size> ws = { 0,-1, 0, 0, 0, 0, 1,-1,-1,-1,-1, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 1, 1, 1, 1,-2, 0, 0, 0, 0, 2,-2,-2,-2,-2,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2,-2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,-2,-2,-2,-2, 0, 0, 0, 0, 2, 2, 2, 2,-2,-2,-2,-2,-2,-2,-2,-2,-1,-1,-1,-1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2};

	auto sq_dist = std::numeric_limits<float>::max();
	auto cell_id = 0;
	glm::vec3 point;

	for (auto loop = 0; loop < size; ++loop)
	{
		const auto shift  = glm::ivec3(us[loop], vs[loop], ws[loop]);
		const auto hash   = my_hash(quantized + shift);
		const auto offset = glm::vec3(OffsetX(hash, jitter), OffsetY(hash, jitter), OffsetZ(hash, jitter));
		const auto sample = origin + offset + glm::vec3(shift);
		const auto tmp    = glm::dot(source - sample, source - sample);
		if (sq_dist > tmp)
		{
			sq_dist = tmp;
			cell_id = hash;
			point   = sample;
		}
	}

	return std::make_tuple(cell_id, sq_dist, point);
}

// octant optimization
std::tuple<int32_t, float, glm::vec3> Voronoi::Evaluate3DOpt(const glm::vec3& source, int32_t (*my_hash)(const glm::ivec3& p), const float jitter)
{
	assert(0.0f <= jitter && jitter <= 1.0f);

	const auto origin    = glm::vec3(std::floor(source.x), std::floor(source.y), std::floor(source.z));
	const auto quantized = glm::ivec3(int32_t(origin.x), int32_t(origin.y), int32_t(origin.z));

	// early termination
	if (0.0f == jitter)
	{
		const auto hash   = my_hash(quantized);
		const auto sample = origin + 0.5f;
		return std::make_tuple(hash, glm::dot(source - sample, source - sample), sample);
	}

	// octant selection
	auto octant = 0;
	if (0.5f > source.x - origin.x) { octant += 1; }
	if (0.5f > source.y - origin.y) { octant += 2; }
	if (0.5f > source.z - origin.z) { octant += 4; }

	auto sq_dist = std::numeric_limits<float>::max();
	auto cell_id = 0;
	glm::vec3 point;

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
					const auto shift  = glm::ivec3(us3[octant][loop], vs3[octant][loop], ws3[octant][loop]);
					const auto hash   = my_hash(quantized + shift);
					const auto offset = glm::vec3(OffsetX(hash, jitter), OffsetY(hash, jitter), OffsetZ(hash, jitter));
					const auto sample = origin + offset + glm::vec3(shift);
					const auto tmp    = glm::dot(source - sample, source - sample);
					if (sq_dist > tmp)
					{
						sq_dist = tmp;
						cell_id = hash;
						point   = sample;
					}
				}
			}
			else
			{
				break;
			}
		}

		return std::make_tuple(cell_id, sq_dist, point);
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
				const auto shift  = glm::ivec3(us3[octant][loop], vs3[octant][loop], ws3[octant][loop]);
				const auto hash   = my_hash(quantized + shift);
				const auto offset = glm::vec3(OffsetX(hash, jitter), OffsetY(hash, jitter), OffsetZ(hash, jitter));
				const auto sample = origin + offset + glm::vec3(shift);
				const auto tmp    = glm::dot(source - sample, source - sample);
				if (sq_dist > tmp)
				{
					sq_dist = tmp;
					cell_id = hash;
					point   = sample;
				}
			}
		}
		else
		{
			break;
		}
	}

	return std::make_tuple(cell_id, sq_dist, point);
}

// octant optimization with cache
std::tuple<int32_t, float, glm::vec3> Voronoi::Evaluate3DCache(const glm::vec3& source, int32_t(*my_hash)(const glm::ivec3& p), const float jitter)
{
	assert(0.0f <= jitter && jitter <= 1.0f);

	const auto origin    = glm::vec3(std::floor(source.x), std::floor(source.y), std::floor(source.z));
	const auto quantized = glm::ivec3(int32_t(origin.x), int32_t(origin.y), int32_t(origin.z));
	const auto self      = my_hash(quantized);

	// early termination
	if (0.0f == jitter)
	{
		const auto sample = origin + 0.5f;
		return std::make_tuple(self, glm::dot(source - sample, source - sample), sample);
	}

	// octant selection
	auto octant = 0;
	if (0.5f > source.x - origin.x) { octant += 1; }
	if (0.5f > source.y - origin.y) { octant += 2; }
	if (0.5f > source.z - origin.z) { octant += 4; }

	auto sq_dist = std::numeric_limits<float>::max();
	auto cell_id = 0;
	glm::vec3 point;

	// fast enough
	if (0.5f >= jitter)
	{
		// 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19
		// 0, 1, 1, 1, 2, 2, 2, 3, 9, 9, 9,10,10,10,10,10,10,11,11,11
		const std::array<int32_t, 7> ranges = { 0, 1, 2, 3, 9, 10, 11 };
		const std::array<int32_t, 8> slices = { 0, 1, 4, 7, 8, 11, 17, 20 };

		if (cache3d_cell_id == self && cache3d_octant == octant)
		{
			// cache hit
			auto counter = 0;
			for (auto dist = 0; dist < 7; ++dist)
			{
				if (ranges[dist] < sq_dist * 16)
				{
					for (auto loop = slices[dist]; loop < slices[dist + 1]; ++loop)
					{
						if (counter < cache3d_counter)
						{
							// reuse cache
							const auto sample = cache3d_samples[counter];
							const auto tmp    = glm::dot(source - sample, source - sample);
							if (sq_dist > tmp)
							{
								sq_dist = tmp;
								cell_id = cache3d_cell_ids[counter];
								point   = sample;
							}
						}
						else
						{
							// fall back
							const auto shift  = glm::ivec3(us3[octant][loop], vs3[octant][loop], ws3[octant][loop]);
							const auto hash   = my_hash(quantized + shift);
							const auto offset = glm::vec3(OffsetX(hash, jitter), OffsetY(hash, jitter), OffsetZ(hash, jitter));
							const auto sample = origin + offset + glm::vec3(shift);
							const auto tmp    = glm::dot(source - sample, source - sample);

							if (sq_dist > tmp)
							{
								sq_dist = tmp;
								cell_id = hash;
								point   = sample;
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
			cache3d_cell_id = self;
			cache3d_octant  = octant;
			cache3d_counter = 0;

			for (auto dist = 0; dist < 7; ++dist)
			{
				if (ranges[dist] < sq_dist * 16)
				{
					for (auto loop = slices[dist]; loop < slices[dist + 1]; ++loop)
					{
						const auto shift  = glm::ivec3(us3[octant][loop], vs3[octant][loop], ws3[octant][loop]);
						const auto hash   = my_hash(quantized + shift);
						const auto offset = glm::vec3(OffsetX(hash, jitter), OffsetY(hash, jitter), OffsetZ(hash, jitter));
						const auto sample = origin + offset + glm::vec3(shift);
						const auto tmp    = glm::dot(source - sample, source - sample);

						// update cache
						cache3d_samples [cache3d_counter] = sample;
						cache3d_cell_ids[cache3d_counter] = hash;
						++cache3d_counter;

						if (sq_dist > tmp)
						{
							sq_dist = tmp;
							cell_id = hash;
							point   = sample;
						}
					}
				}
				else
				{
					break;
				}
			}
		}

		return std::make_tuple(cell_id, sq_dist, point);
	}

	// cell id and possible minimum squared distance
	// 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38
	// 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
	const std::array<int32_t, 5> ranges = { 0,  1,  2,  3,  4 };
	const std::array<int32_t, 6> slices = { 0,  8, 20, 26, 27, 39 };

	if (cache3d_cell_id == self && cache3d_octant == octant)
	{
		// cache hit
		auto counter = 0;
		for (auto dist = 0; dist < 5; ++dist)
		{
			if (ranges[dist] < sq_dist * 4)
			{
				for (auto loop = slices[dist]; loop < slices[dist + 1]; ++loop)
				{
					if (counter < cache3d_counter)
					{
						// reuse cache
						const auto sample = cache3d_samples[counter];
						const auto tmp    = glm::dot(source - sample, source - sample);
						if (sq_dist > tmp)
						{
							sq_dist = tmp;
							cell_id = cache3d_cell_ids[counter];
							point   = sample;
						}
					}
					else
					{
						// fall back
						const auto shift  = glm::ivec3(us3[octant][loop], vs3[octant][loop], ws3[octant][loop]);
						const auto hash   = my_hash(quantized + shift);
						const auto offset = glm::vec3(OffsetX(hash, jitter), OffsetY(hash, jitter), OffsetZ(hash, jitter));
						const auto sample = origin + offset + glm::vec3(shift);
						const auto tmp    = glm::dot(source - sample, source - sample);
						
						if (sq_dist > tmp)
						{
							sq_dist = tmp;
							cell_id = hash;
							point   = sample;
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
		cache3d_cell_id = self;
		cache3d_octant  = octant;
		cache3d_counter = 0;

		for (auto dist = 0; dist < 5; ++dist)
		{
			if (ranges[dist] < sq_dist * 4)
			{
				for (auto loop = slices[dist]; loop < slices[dist + 1]; ++loop)
				{
					const auto shift  = glm::ivec3(us3[octant][loop], vs3[octant][loop], ws3[octant][loop]);
					const auto hash   = my_hash(quantized + shift);
					const auto offset = glm::vec3(OffsetX(hash, jitter), OffsetY(hash, jitter), OffsetZ(hash, jitter));
					const auto sample = origin + offset + glm::vec3(shift);
					const auto tmp    = glm::dot(source - sample, source - sample);

					// update cache
					cache3d_samples [cache3d_counter] = sample;
					cache3d_cell_ids[cache3d_counter] = hash;
					++cache3d_counter;

					if (sq_dist > tmp)
					{
						sq_dist = tmp;
						cell_id = hash;
						point   = sample;
					}
				}
			}
			else
			{
				break;
			}
		}
	}

	return std::make_tuple(cell_id, sq_dist, point);
}

// naive implementation
std::tuple<int32_t, float, glm::vec4> Voronoi::Evaluate4DRef(const glm::vec4& source, int32_t (*my_hash)(const glm::ivec4& p), const float jitter)
{
	assert(0.0f <= jitter && jitter <= 1.0f);

	const auto origin    = glm::vec4(std::floor(source.x), std::floor(source.y), std::floor(source.z), std::floor(source.w));
	const auto quantized = glm::ivec4(int32_t(origin.x), int32_t(origin.y), int32_t(origin.z), int32_t(origin.w));

	auto sq_dist = std::numeric_limits<float>::max();
	auto cell_id = 0;
	glm::vec4 point;

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
					const auto shift  = glm::ivec4(u, v, w, t);
					const auto hash   = my_hash(quantized + shift);
					const auto offset = glm::vec4(OffsetX(hash, jitter), OffsetY(hash, jitter), OffsetZ(hash, jitter) ,OffsetT(hash, jitter));
					const auto sample = origin + offset + glm::vec4(shift);
					const auto tmp    = glm::dot(source - sample, source - sample);
					if (sq_dist > tmp)
					{
						sq_dist = tmp;
						cell_id = hash;
						point   = sample;

					}
				}
			}
		}
	}

	return std::make_tuple(cell_id, sq_dist, point);
}

// hextant optimization
std::tuple<int32_t, float, glm::vec4> Voronoi::Evaluate4DOpt(const glm::vec4& source, int32_t(*my_hash)(const glm::ivec4& p), const float jitter)
{
	assert(0.0f <= jitter && jitter <= 1.0f);

	const auto origin    = glm::vec4(std::floor(source.x), std::floor(source.y), std::floor(source.z), std::floor(source.w));
	const auto quantized = glm::ivec4(int32_t(origin.x), int32_t(origin.y), int32_t(origin.z), int32_t(origin.w));

	// early termination
	if (0.0f == jitter)
	{
		const auto hash   = my_hash(quantized);
		const auto sample = origin + 0.5f;
		return std::make_tuple(hash, glm::dot(source - sample, source - sample), sample);
	}

	// hextant selection
	auto hextant = 0;
	if (0.5f > source.x - origin.x) { hextant += 1; }
	if (0.5f > source.y - origin.y) { hextant += 2; }
	if (0.5f > source.z - origin.z) { hextant += 4; }
	if (0.5f > source.t - origin.t) { hextant += 8; }

	auto sq_dist = std::numeric_limits<float>::max();
	auto cell_id = 0;
	glm::vec4 point;

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
					const auto shift  = glm::ivec4(us4[hextant][loop], vs4[hextant][loop], ws4[hextant][loop], ts4[hextant][loop]);
					const auto hash   = my_hash(quantized + shift);
					const auto offset = glm::vec4(OffsetX(hash, jitter), OffsetY(hash, jitter), OffsetZ(hash, jitter) ,OffsetT(hash, jitter));
					const auto sample = origin + offset + glm::vec4(shift);
					const auto tmp    = glm::dot(source - sample, source - sample);
					if (sq_dist > tmp)
					{
						sq_dist = tmp;
						cell_id = hash;
						point   = sample;
					}
				}
			}
			else
			{
				break;
			}
		}

		return std::make_tuple(cell_id, sq_dist, point);
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
				const auto shift  = glm::ivec4(us4[hextant][loop], vs4[hextant][loop], ws4[hextant][loop], ts4[hextant][loop]);
				const auto hash   = my_hash(quantized + shift);
				const auto offset = glm::vec4(OffsetX(hash, jitter), OffsetY(hash, jitter), OffsetZ(hash, jitter) ,OffsetT(hash, jitter));
				const auto sample = origin + offset + glm::vec4(shift);
				const auto tmp    = glm::dot(source - sample, source - sample);
				if (sq_dist > tmp)
				{
					sq_dist = tmp;
					cell_id = hash;
					point   = sample;
				}
			}
		}
		else
		{
			break;
		}
	}

	return std::make_tuple(cell_id, sq_dist, point);
}

// hextant optimization with cache
std::tuple<int32_t, float, glm::vec4> Voronoi::Evaluate4DCache(const glm::vec4& source, int32_t(*my_hash)(const glm::ivec4& p), const float jitter)
{
	assert(0.0f <= jitter && jitter <= 1.0f);

	const auto origin    = glm::vec4(std::floor(source.x), std::floor(source.y), std::floor(source.z), std::floor(source.w));
	const auto quantized = glm::ivec4(int32_t(origin.x), int32_t(origin.y), int32_t(origin.z), int32_t(origin.w));
	const auto self      = my_hash(quantized);

	// early termination
	if (0.0f == jitter)
	{
		const auto sample = origin + 0.5f;
		return std::make_tuple(self, glm::dot(source - sample, source - sample), sample);
	}

	// hextant selection
	auto hextant = 0;
	if (0.5f > source.x - origin.x) { hextant += 1; }
	if (0.5f > source.y - origin.y) { hextant += 2; }
	if (0.5f > source.z - origin.z) { hextant += 4; }
	if (0.5f > source.t - origin.t) { hextant += 8; }

	auto sq_dist = std::numeric_limits<float>::max();
	auto cell_id = 0;
	glm::vec4 point;

	if (0.5f >= jitter)
	{
		// cell id and possible minimum squared distance
		//  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84
		//  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10
		const std::array<int32_t, 6> ranges = { 0, 1,  2,  3,  9, 10 };
		const std::array<int32_t, 7> slices = { 0, 5, 20, 25, 40, 55, 85 };

		if (cache4d_cell_id == self && cache4d_hextant == hextant)
		{
			// cache hit
			auto counter = 0;
			for (auto dist = 0; dist < 6; ++dist)
			{
				if (ranges[dist] < sq_dist * 16)
				{
					for (auto loop = slices[dist]; loop < slices[dist + 1]; ++loop)
					{
						if (counter < cache4d_counter)
						{
							// reuse cache
							const auto sample = cache4d_samples[counter];
							const auto tmp    = glm::dot(source - sample, source - sample);
							if (sq_dist > tmp)
							{
								sq_dist = tmp;
								cell_id = cache4d_cell_ids[counter];
								point   = sample;
							}
						}
						else
						{
							// fall back
							const auto shift  = glm::ivec4(us4[hextant][loop], vs4[hextant][loop], ws4[hextant][loop], ts4[hextant][loop]);
							const auto hash   = my_hash(quantized + shift);
							const auto offset = glm::vec4(OffsetX(hash, jitter), OffsetY(hash, jitter), OffsetZ(hash, jitter), OffsetT(hash, jitter));
							const auto sample = origin + offset + glm::vec4(shift);
							const auto tmp    = glm::dot(source - sample, source - sample);

							if (sq_dist > tmp)
							{
								sq_dist = tmp;
								cell_id = hash;
								point   = sample;
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
			cache4d_cell_id = self;
			cache4d_hextant = hextant;
			cache4d_counter = 0;

			for (auto dist = 0; dist < 6; ++dist)
			{
				if (ranges[dist] < sq_dist * 16)
				{
					for (auto loop = slices[dist]; loop < slices[dist + 1]; ++loop)
					{
						const auto shift  = glm::ivec4(us4[hextant][loop], vs4[hextant][loop], ws4[hextant][loop], ts4[hextant][loop]);
						const auto hash   = my_hash(quantized + shift);
						const auto offset = glm::vec4(OffsetX(hash, jitter), OffsetY(hash, jitter), OffsetZ(hash, jitter), OffsetT(hash, jitter));
						const auto sample = origin + offset + glm::vec4(shift);
						const auto tmp    = glm::dot(source - sample, source - sample);

						// update cache
						cache4d_samples [cache4d_counter] = sample;
						cache4d_cell_ids[cache4d_counter] = hash;
						++cache4d_counter;

						if (sq_dist > tmp)
						{
							sq_dist = tmp;
							cell_id = hash;
							point   = sample;
						}
					}
				}
				else
				{
					break;
				}
			}
		}

		return std::make_tuple(cell_id, sq_dist, point);
	}

	// cell id and possible minimum squared distance
	// 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194
	// 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  3,  3,  3,  3,  3,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4
	const std::array<int32_t, 5> ranges = { 0,   1,   2,   3,   4 };
	const std::array<int32_t, 6> slices = { 0,  40, 100, 130, 135, 195 };

	if (cache4d_cell_id == self && cache4d_hextant == hextant)
	{
		// cache hit
		auto counter = 0;
		for (auto dist = 0; dist < 5; ++dist)
		{
			if (ranges[dist] < sq_dist * 4)
			{
				for (auto loop = slices[dist]; loop < slices[dist + 1]; ++loop)
				{
					if (counter < cache4d_counter)
					{
						// reuse cache
						const auto sample = cache4d_samples[counter];
						const auto tmp    = glm::dot(source - sample, source - sample);
						if (sq_dist > tmp)
						{
							sq_dist = tmp;
							cell_id = cache4d_cell_ids[counter];
							point   = sample;
						}
					}
					else
					{
						// fall back
						const auto shift  = glm::ivec4(us4[hextant][loop], vs4[hextant][loop], ws4[hextant][loop], ts4[hextant][loop]);
						const auto hash   = my_hash(quantized + shift);
						const auto offset = glm::vec4(OffsetX(hash, jitter), OffsetY(hash, jitter), OffsetZ(hash, jitter), OffsetT(hash, jitter));
						const auto sample = origin + offset + glm::vec4(shift);
						const auto tmp    = glm::dot(source - sample, source - sample);

						if (sq_dist > tmp)
						{
							sq_dist = tmp;
							cell_id = hash;
							point   = sample;
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
		cache4d_cell_id = self;
		cache4d_hextant = hextant;
		cache4d_counter = 0;

		for (auto dist = 0; dist < 5; ++dist)
		{
			if (ranges[dist] < sq_dist * 4)
			{
				for (auto loop = slices[dist]; loop < slices[dist + 1]; ++loop)
				{
					const auto shift  = glm::ivec4(us4[hextant][loop], vs4[hextant][loop], ws4[hextant][loop], ts4[hextant][loop]);
					const auto hash   = my_hash(quantized + shift);
					const auto offset = glm::vec4(OffsetX(hash, jitter), OffsetY(hash, jitter), OffsetZ(hash, jitter), OffsetT(hash, jitter));
					const auto sample = origin + offset + glm::vec4(shift);
					const auto tmp    = glm::dot(source - sample, source - sample);

					// update cache
					cache4d_samples [cache4d_counter] = sample;
					cache4d_cell_ids[cache4d_counter] = hash;
					++cache4d_counter;

					if (sq_dist > tmp)
					{
						sq_dist = tmp;
						cell_id = hash;
						point   = sample;
					}
				}
			}
			else
			{
				break;
			}
		}
	}

	return std::make_tuple(cell_id, sq_dist, point);
}

// naive triangle implementation
std::tuple<int32_t, float, glm::vec2> Voronoi::Evaluate2DTri(const glm::vec2& source, int32_t(*my_hash)(const glm::ivec2& p), const float jitter)
{
	assert(0.0f <= jitter && jitter <= 1.0f);

	const auto local = glm::vec2(glm::dot(glm::vec2(1.0f, -one_over_sqrt3), source), glm::dot(glm::vec2(0.0f, 2.0f * one_over_sqrt3), source));

	const auto origin    = glm::vec2(std::floor(local.x), std::floor(local.y));
	const auto quantized = glm::ivec2(int32_t(origin.x), int32_t(origin.y));

	auto sq_dist = std::numeric_limits<float>::max();
	auto cell_id = 0;
	glm::vec2 point;

	// 1 (self) + 14 (neighbours)
	const auto size = 15;
	const std::array<int32_t, size> us[2] =
	{
		{ 0, 0,-1, 1, 0,-1, 1,-1, 1, -1,-2,-2, 2, 1, 0 },
		{ 0, 0,-1, 1, 0,-1, 1,-1, 1,  0,-1,-2, 2, 2, 1 },
	};
	const std::array<int32_t, size> vs[2] =
	{
		{ 0,-1, 0, 0, 1,-1,-1, 1, 1, 2, 1, 0,-1,-2,-2 },
		{ 0,-1, 0, 0, 1,-1,-1, 1, 1, 2, 2, 1, 0,-1,-2 },
	};
	const std::array<int32_t, size> flags[2] =
	{
		{ 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 3, 2, 1, 3, 2 },
		{ 3, 3, 3, 3, 3, 2, 3, 3, 3, 1, 3, 2, 1, 3, 2 },
	};

	// "A Low-Distortion Map Between Triangle and Square" by Eric Heitz
	// maps a unit - square point (x, y) to a unit - triangle point
	auto triangle = [&](float& x, float& y)
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

	const auto tmp   = local - origin;
	const auto which = (1.0f > tmp.x + tmp.y) ? 0 : 1;
	for (auto loop = 0; loop < size; ++loop)
	{
		const auto shift = glm::ivec2(us[which][loop], vs[which][loop]);

		// lower triangle
		if (1 & flags[which][loop])
		{
			const auto hash = my_hash(quantized + shift + 0);
			auto randomX = OffsetX(hash);
			auto randomY = OffsetY(hash);
			triangle(randomX, randomY);
			const auto offset = glm::vec2(randomX, randomY) * jitter + 1.0f / 3.0f;
			const auto sample = origin + offset + glm::vec2(shift);
			const auto global = glm::vec2(glm::dot(glm::vec2(1.0f, 0.5f), sample), glm::dot(glm::vec2(0.0f, std::sqrt(3) * 0.5f), sample));
			const auto tmp    = glm::dot(source - global, source - global);
			if (sq_dist > tmp)
			{
				sq_dist = tmp;
				cell_id = hash;
				point   = sample;
			}
		}

		// upper triangle
		if (2 & flags[which][loop])
		{
			const auto hash = my_hash(quantized + shift + PrimeW);
			auto randomX = OffsetX(hash);
			auto randomY = OffsetY(hash);
			triangle(randomX, randomY);
			const auto offset = (glm::vec2(0.5f, 0.5f) - glm::vec2(randomX, randomY)) * jitter + 2.0f / 3.0f;
			const auto sample = origin + offset + glm::vec2(shift);
			const auto global = glm::vec2(glm::dot(glm::vec2(1.0f, 0.5f), sample), glm::dot(glm::vec2(0.0f, std::sqrt(3) * 0.5f), sample));
			const auto tmp    = glm::dot(source - global, source - global);
			if (sq_dist > tmp)
			{
				sq_dist = tmp;
				cell_id = hash;
				point   = sample;
			}
		}
	}

	return std::make_tuple(cell_id, sq_dist, point);
}

// naive honeycomb implementation
std::tuple<int32_t, float, glm::vec2> Voronoi::Evaluate2DHex(const glm::vec2& source, int32_t(*my_hash)(const glm::ivec2& p), const float jitter)
{
	assert(0.0f <= jitter && jitter <= 1.0f);

	static const std::array<glm::vec2, 3> basis =
	{
		glm::vec2(-sqrt3 * 0.5f, -0.5f),
		glm::vec2( sqrt3 * 0.5f, -0.5f),
		glm::vec2(0, 1)
	};

	static const std::array<glm::vec2, 3> ortho =
	{
		glm::vec2(   -one_over_sqrt3, -1),
		glm::vec2(    one_over_sqrt3, -1),
		glm::vec2(2 * one_over_sqrt3,  0)
	};

	const std::array<float, 3> dots =
	{
		glm::dot(ortho[0], source),
		glm::dot(ortho[1], source),
		glm::dot(ortho[2], source)
	};

	const std::array<float, 6> grids =
	{
		std::floor( dots[0]), std::floor( dots[1]), std::floor( dots[2]),
		std::floor(-dots[0]), std::floor(-dots[1]), std::floor(-dots[2])
	};

	const int32_t int_grids[6] =
	{
		int32_t(grids[0]), int32_t(grids[1]), int32_t(grids[2]),
		int32_t(grids[3]), int32_t(grids[4]), int32_t(grids[5])
	};

	glm::ivec2 quantized;
	glm::vec2  origin(0, 0);
	if ((int_grids[0] + int_grids[1]) % 3 == 0)
	{
		origin    = basis[0] * grids[0] + basis[1] * grids[1];
		quantized = { int_grids[0], int_grids[1] };
	}
	else if ((int_grids[2] + int_grids[3]) % 3 == 0)
	{
		origin    = basis[1] * grids[2] + basis[2] * grids[3];
		quantized = { -int_grids[3], int_grids[2] - int_grids[3] };
	}
	else if ((int_grids[4] + int_grids[5]) % 3 == 0)
	{
		origin    = basis[2] * grids[4] + basis[0] * grids[5];
		quantized = { int_grids[5] - int_grids[4], -int_grids[4] };
	}

	const std::array<int32_t, 4>  slices = { 0, 7, 13, 19 };
	const std::array<  float, 3>  ranges = { 0, 1, 3 };
	const std::array<int32_t, 19> us     = { 0, 1, -1, -2, -1,  1, 2, 3, 0, -3, -3,  0, 3, 2, -2, -4, -2,  2, 4 };
	const std::array<int32_t, 19> vs     = { 0, 2,  1, -1, -2, -1, 1, 3, 3,  0, -3, -3, 0, 4,  2, -2, -4, -2, 2 };

	auto hexagon = [&](const int32_t hash)
	{
		const auto axis = hash % 3;
		return basis[axis] * OffsetX(hash) + basis[(axis+1) % 3] * OffsetY(hash);
	};

	auto sq_dist = std::numeric_limits<float>::max();
	auto cell_id = 0;
	glm::vec2 point;

	for (auto dist = 0; dist < 3; ++dist)
	{
		if (ranges[dist] < sq_dist)
		{
			for (auto loop = slices[dist]; loop < slices[dist + 1]; ++loop)
			{
				const auto hash   = my_hash(quantized + glm::ivec2(us[loop], vs[loop]));
				const auto offset = hexagon(hash) * jitter;
				const auto sample = origin + offset + basis[0] * float(us[loop]) + basis[1] * float(vs[loop]);
				const auto tmp    = glm::dot(source - sample, source - sample);
				if (sq_dist > tmp)
				{
					sq_dist = tmp;
					cell_id = hash;
					point   = sample;
				}
			}
		}
		else
		{
			break;
		}
	}

	return std::make_tuple(cell_id, sq_dist, point);
}

const int32_t Voronoi::us2[4][13] =
{
	{ 0, 1, 0, 1, 0,-1, 1,-1,-1, 2, 0, 2, 1},
	{ 0,-1, 0,-1, 0, 1,-1, 1, 1,-2, 0,-2,-1},
	{ 0, 0, 1, 1,-1, 0,-1, 1,-1, 0, 2, 1, 2},
	{ 0, 0,-1,-1, 1, 0, 1,-1, 1, 0,-2,-1,-2},
};

const int32_t Voronoi::vs2[4][13] =
{
	{ 0, 0, 1, 1,-1, 0,-1, 1,-1, 0, 2, 1, 2},
	{ 0, 0, 1, 1,-1, 0,-1, 1,-1, 0, 2, 1, 2},
	{ 0,-1, 0,-1, 0, 1,-1, 1, 1,-2, 0,-2,-1},
	{ 0,-1, 0,-1, 0, 1,-1, 1, 1,-2, 0,-2,-1},
};

const int32_t Voronoi::us3[8][39] =
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

const int32_t Voronoi::vs3[8][39] =
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

const int32_t Voronoi::ws3[8][39] =
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

const int32_t Voronoi::us4[16][195] =
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

const int32_t Voronoi::vs4[16][195] =
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

const int32_t Voronoi::ws4[16][195] =
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

const int32_t Voronoi::ts4[16][195] =
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
#endif
}

#endif
