// APACHE LICENSE, VERSION 2.0
// copyright (c) shinji ogaki

#include <algorithm>
#include <chrono>
#include <execution>
#include <iostream>
#include <iomanip>
#include <random>

#define VOROCE_IMPLEMENTATION
#include "../include/voroce.h"

int main()
{
	using voroce::Voronoi;

	const auto golden_ratio = 1.61803398875;
	const auto N = 1000000;
	const auto T = 10000000;

	////////////////////////////////////////////////////////////////
	// 2D
	////////////////////////////////////////////////////////////////
	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto v = Voronoi::Evaluate2DRef(glm::vec2(x * 2048 - 1024, y * 2048 - 1024), Voronoi::Hash2DLowQuality, 0.5);
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "2d ref half: " << elapsed << std::endl;
	}

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto v = Voronoi::Evaluate2DOpt(glm::vec2(x * 2048 - 1024, y * 2048 - 1024), Voronoi::Hash2DLowQuality, 0.5f);
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "2d opt half: " << elapsed << std::endl;
	}

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto v = Voronoi::Evaluate2DRef(glm::vec2(x * 2048 - 1024, y * 2048 - 1024), Voronoi::Hash2DLowQuality);
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "2d ref full: " << elapsed << std::endl;
	}

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto v = Voronoi::Evaluate2DOpt(glm::vec2(x * 2048 - 1024, y * 2048 - 1024), Voronoi::Hash2DLowQuality);
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "2d opt full: " << elapsed << std::endl;
	}

	{
		// bug or numerical error
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto ref = Voronoi::Evaluate2DRef(glm::vec2(x * 2048 - 1024, y * 2048 - 1024), Voronoi::Hash2DLowQuality, 0.5f);
			const auto opt = Voronoi::Evaluate2DOpt(glm::vec2(x * 2048 - 1024, y * 2048 - 1024), Voronoi::Hash2DLowQuality, 0.5f);
			if (ref != opt)
			{
				std::cout << "NG: (" << x << "," << y << ") ref:" << get<0>(ref) << " opt:" << get<0>(opt) << std::endl;
			}
		}
	}

	{
		// bug or numerical error
		std::mt19937 mt(15);
		std::uniform_real_distribution<double> uni(0.0f, 1.0f);
		for (auto L = 0; L < 256; ++L)
		{
			std::cout << "2d: " << L << std::endl;
			std::vector<std::tuple<double, double>> RND(T);
			for (auto& r : RND)
			{
				r = { uni(mt), uni(mt) };
			}
			std::for_each(std::execution::par, RND.begin(), RND.end(), [&](std::tuple<double, double>& tuple)
				{
					const auto [x, y] = tuple;
					const auto ref = Voronoi::Evaluate2DRef(glm::vec2(x * 2048 - 1024, y * 2048 - 1024), Voronoi::Hash2DLowQuality, 0.5f);
					const auto opt = Voronoi::Evaluate2DOpt(glm::vec2(x * 2048 - 1024, y * 2048 - 1024), Voronoi::Hash2DLowQuality, 0.5f);
					if (get<1>(ref) != get<1>(opt))
					{
						std::cout << "2d: ref:" << get<0>(ref) << " opt:" << get<0>(opt) << " ref:" << std::setprecision(16) << get<1>(ref) << " opt:" << get<1>(opt) << std::endl;
					}
				});
		}
	}

	{
		// bug or numerical error
		std::mt19937 mt(15);
		std::uniform_real_distribution<double> uni(0.0f, 1.0f);
		for (auto L = 0; L < 256; ++L)
		{
			std::cout << "2d: " << L << std::endl;
			std::vector<std::tuple<double, double>> RND(T);
			for (auto& r : RND)
			{
				r = { uni(mt), uni(mt) };
			}
			std::for_each(std::execution::par, RND.begin(), RND.end(), [&](std::tuple<double, double>& tuple)
				{
					const auto [x, y] = tuple;
					const auto ref = Voronoi::Evaluate2DRef(glm::vec2(x * 2048 - 1024, y * 2048 - 1024), Voronoi::Hash2DLowQuality);
					const auto opt = Voronoi::Evaluate2DOpt(glm::vec2(x * 2048 - 1024, y * 2048 - 1024), Voronoi::Hash2DLowQuality);
					if (get<1>(ref) != get<1>(opt))
					{
						std::cout << "2d: ref:" << get<0>(ref) << " opt:" << get<0>(opt) << " ref:" << std::setprecision(16) << get<1>(ref) << " opt:" << get<1>(opt) << std::endl;
					}
				});
		}
	}

	////////////////////////////////////////////////////////////////
	// 3D
	////////////////////////////////////////////////////////////////
	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto z = y;
			const auto v = Voronoi::Evaluate3DRef(glm::vec3(x * 2048 - 1024, y * 2048 - 1024, z * 2048 - 1024), Voronoi::Hash3DLowQuality, 0.5f);
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "3d ref half: " << elapsed << std::endl;
	}

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto z = y;
			const auto v = Voronoi::Evaluate3DOpt(glm::vec3(x * 2048 - 1024, y * 2048 - 1024, z * 2048 - 1024), Voronoi::Hash3DLowQuality, 0.5f);
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "3d opt half: " << elapsed << std::endl;
	}

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto z = y;
			const auto v = Voronoi::Evaluate3DRef(glm::vec3(x * 2048 - 1024, y * 2048 - 1024, z * 2048 - 1024), Voronoi::Hash3DLowQuality);
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "3d ref full: " << elapsed << std::endl;
	}

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto z = y;
			const auto v = Voronoi::Evaluate3DOpt(glm::vec3(x * 2048 - 1024, y * 2048 - 1024, z * 2048 - 1024), Voronoi::Hash3DLowQuality);
		}
		const auto end = std::chrono::system_clock::now();
		const auto  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "3d opt full: " << elapsed << std::endl;
	}

	{
		// bug or numerical error
		std::mt19937 mt(15);
		std::uniform_real_distribution<double> uni(0.0f, 1.0f);
		for (auto L = 0; L < 256; ++L)
		{
			std::cout << "3d: " << L << std::endl;
			std::vector<std::tuple<double, double, double>> RND(T);
			for (auto& r : RND)
			{
				r = { uni(mt), uni(mt), uni(mt) };
			}
			std::for_each(std::execution::par, RND.begin(), RND.end(), [&](std::tuple<double, double, double>& tuple)
				{
					const auto [x, y, z] = tuple;
					const auto ref = Voronoi::Evaluate3DRef(glm::vec3(x * 2048 - 1024, y * 2048 - 1024, z * 2048 - 1024), Voronoi::Hash3DLowQuality, 0.5f);
					const auto opt = Voronoi::Evaluate3DOpt(glm::vec3(x * 2048 - 1024, y * 2048 - 1024, z * 2048 - 1024), Voronoi::Hash3DLowQuality, 0.5f);
					if (get<1>(ref) != get<1>(opt))
					{
						std::cout << "3d: ref:" << get<0>(ref) << " opt:" << get<0>(opt) << " ref:" << std::setprecision(16) << get<1>(ref) << " opt:" << get<1>(opt) << std::endl;
					}
				});
		}
	}

	{
		// bug or numerical error
		std::mt19937 mt(15);
		std::uniform_real_distribution<double> uni(0.0f, 1.0f);
		for (auto L = 0; L < 256; ++L)
		{
			std::cout << "3d: " << L << std::endl;
			std::vector<std::tuple<double, double, double>> RND(T);
			for (auto& r : RND)
			{
				r = { uni(mt), uni(mt), uni(mt) };
			}
			std::for_each(std::execution::par, RND.begin(), RND.end(), [&](std::tuple<double, double, double>& tuple)
				{
					const auto [x, y, z] = tuple;
					const auto ref = Voronoi::Evaluate3DRef(glm::vec3(x * 2048 - 1024, y * 2048 - 1024, z * 2048 - 1024), Voronoi::Hash3DLowQuality);
					const auto opt = Voronoi::Evaluate3DOpt(glm::vec3(x * 2048 - 1024, y * 2048 - 1024, z * 2048 - 1024), Voronoi::Hash3DLowQuality);
					if (get<1>(ref) != get<1>(opt))
					{
						std::cout << "3d: ref:" << get<0>(ref) << " opt:" << get<0>(opt) << " ref:" << std::setprecision(16) << get<1>(ref) << " opt:" << get<1>(opt) << std::endl;
					}
				});
		}
	}

	////////////////////////////////////////////////////////////////
	// 4D
	////////////////////////////////////////////////////////////////
	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto z = y;
			const auto t = x;
			const auto v = Voronoi::Evaluate4DRef(glm::vec4(x * 2048 - 1024, y * 2048 - 1024, z * 2048 - 1024, t * 2048 - 1024), Voronoi::Hash4DLowQuality, 0.5f);
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "4d ref half: " << elapsed << std::endl;
	}

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto z = y;
			const auto t = x;
			const auto v = Voronoi::Evaluate4DOpt(glm::vec4(x * 2048 - 1024, y * 2048 - 1024, z * 2048 - 1024, t * 2048 - 1024), Voronoi::Hash4DLowQuality, 0.5f);
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "4d opt half: " << elapsed << std::endl;
	}

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto z = y;
			const auto t = x;
			const auto v = Voronoi::Evaluate4DRef(glm::vec4(x * 2048 - 1024, y * 2048 - 1024, z * 2048 - 1024, t * 2048 - 1024), Voronoi::Hash4DLowQuality);
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "4d ref full: " << elapsed << std::endl;
	}

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto z = y;
			const auto t = x;
			const auto v = Voronoi::Evaluate4DOpt(glm::vec4(x * 2048 - 1024, y * 2048 - 1024, z * 2048 - 1024, t * 2048 - 1024), Voronoi::Hash4DLowQuality);
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "4d opt full: " << elapsed << std::endl;
	}

	{
		// bug or numerical error
		std::mt19937 mt(15);
		std::uniform_real_distribution<double> uni(0.0f, 1.0f);
		for (auto L = 0; L < 256; ++L)
		{
			std::cout << "4d: " << L << std::endl;
			std::vector<std::tuple<double, double, double, double>> RND(T);
			for (auto& r : RND)
			{
				r = { uni(mt), uni(mt), uni(mt), uni(mt) };
			}
			std::for_each(std::execution::par, RND.begin(), RND.end(), [&](std::tuple<double, double, double, double>& tuple)
				{
					const auto [x, y, z, t] = tuple;
					const auto ref = Voronoi::Evaluate4DRef(glm::vec4(x * 2048 - 1024, y * 2048 - 1024, z * 2048 - 1024, t * 2048 - 1024), Voronoi::Hash4DLowQuality, 0.5f);
					const auto opt = Voronoi::Evaluate4DOpt(glm::vec4(x * 2048 - 1024, y * 2048 - 1024, z * 2048 - 1024, t * 2048 - 1024), Voronoi::Hash4DLowQuality, 0.5f);
					if (get<1>(ref) != get<1>(opt))
					{
						std::cout << "4d: ref:" << get<0>(ref) << " opt:" << get<0>(opt) << " ref:" << std::setprecision(16) << get<1>(ref) << " opt:" << get<1>(opt) << std::endl;
					}
				});
		}
	}

	{
		// bug or numerical error
		std::mt19937 mt(15);
		std::uniform_real_distribution<double> uni(0.0f, 1.0f);
		for (auto L = 0; L < 256; ++L)
		{
			std::cout << "4d: " << L << std::endl;
			std::vector<std::tuple<double, double, double, double>> RND(T);
			for (auto& r : RND)
			{
				r = { uni(mt), uni(mt), uni(mt), uni(mt) };
			}
			std::for_each(std::execution::par, RND.begin(), RND.end(), [&](std::tuple<double, double, double, double>& tuple)
				{
					const auto [x, y, z, t] = tuple;
					const auto ref = Voronoi::Evaluate4DRef(glm::vec4(x * 2048 - 1024, y * 2048 - 1024, z * 2048 - 1024, t * 2048 - 1024), Voronoi::Hash4DLowQuality);
					const auto opt = Voronoi::Evaluate4DOpt(glm::vec4(x * 2048 - 1024, y * 2048 - 1024, z * 2048 - 1024, t * 2048 - 1024), Voronoi::Hash4DLowQuality);
					if (get<1>(ref) != get<1>(opt))
					{
						std::cout << "4d: ref:" << get<0>(ref) << " opt:" << get<0>(opt) << " ref:" << std::setprecision(16) << get<1>(ref) << " opt:" << get<1>(opt) << std::endl;
					}
				});
		}
	}
}
