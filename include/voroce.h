// APACHE LICENSE, VERSION 2.0
// copyright (c) shinji ogaki

#include "glm/glm.hpp"

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

		// reference implementation
		static std::tuple<int32_t, float, glm::vec2> Evaluate2DRef(const glm::vec2& source, int32_t(*my_hash)(const glm::ivec2& p), const float jitter = 1.0f);
		static std::tuple<int32_t, float, glm::vec3> Evaluate3DRef(const glm::vec3& source, int32_t(*my_hash)(const glm::ivec3& p), const float jitter = 1.0f);
		static std::tuple<int32_t, float, glm::vec4> Evaluate4DRef(const glm::vec4& source, int32_t(*my_hash)(const glm::ivec4& p), const float jitter = 1.0f);

		// optimized (for bounced rays)
		static std::tuple<int32_t, float, glm::vec2> Evaluate2DOpt(const glm::vec2& source, int32_t(*my_hash)(const glm::ivec2& p), const float jitter = 1.0f);
		static std::tuple<int32_t, float, glm::vec3> Evaluate3DOpt(const glm::vec3& source, int32_t(*my_hash)(const glm::ivec3& p), const float jitter = 1.0f);
		static std::tuple<int32_t, float, glm::vec4> Evaluate4DOpt(const glm::vec4& source, int32_t(*my_hash)(const glm::ivec4& p), const float jitter = 1.0f);

		// cached (for primary rays)
		std::tuple<int32_t, float, glm::vec2> Evaluate2DCache(const glm::vec2& source, int32_t(*my_hash)(const glm::ivec2& p), const float jitter = 1.0f);
		std::tuple<int32_t, float, glm::vec3> Evaluate3DCache(const glm::vec3& source, int32_t(*my_hash)(const glm::ivec3& p), const float jitter = 1.0f);
		std::tuple<int32_t, float, glm::vec4> Evaluate4DCache(const glm::vec4& source, int32_t(*my_hash)(const glm::ivec4& p), const float jitter = 1.0f);

		// non rectangular grids
		static std::tuple<int32_t, float, glm::vec2> Evaluate2DTri(const glm::vec2& source, int32_t(*my_hash)(const glm::ivec2& p), const float jitter = 1.0f);
	};
}