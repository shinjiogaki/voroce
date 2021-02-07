// APACHE LICENSE, VERSION 2.0
// copyright (c) shinji ogaki

#include <algorithm>
#include <chrono>
#include <execution>
#include <iostream>
#include <iomanip>
#include <random>

#include "../include/voroce.h"

int main()
{
	const auto golden_ratio = 1.61803398875;
	const auto N = 1000000;
	const auto T = 10000000;

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			volatile auto v = Voroce::Evaluate2DRef(glm::vec2(x * 128 - 64, y * 128 - 64), Voroce::Hash2DLowQuality).first;
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "2d ref: " << elapsed << std::endl;
	}

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			volatile auto v = Voroce::Evaluate2DOpt(glm::vec2(x * 128 - 64, y * 128 - 64), Voroce::Hash2DLowQuality).first;
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "2d opt: " << elapsed << std::endl;
	}

	{
		// bug or numerical error
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto ref = Voroce::Evaluate2DRef(glm::vec2(x * 128 - 64, y * 128 - 64), Voroce::Hash2DLowQuality).first;
			const auto opt = Voroce::Evaluate2DOpt(glm::vec2(x * 128 - 64, y * 128 - 64), Voroce::Hash2DLowQuality).first;
			if (ref != opt)
			{
				std::cout << "NG: (" << x << "," << y << ") ref:" << ref << " opt:" << opt << std::endl;
			}
		}
	}

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto z = y;
			volatile auto v = Voroce::Evaluate3DRef(glm::vec3(x * 128 - 64, y * 128 - 64, z * 128 - 64), Voroce::Hash3DLowQuality).first;
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "3d ref: " << elapsed << std::endl;
	}

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto z = y;
			volatile auto v = Voroce::Evaluate3DOpt(glm::vec3(x * 128 - 64, y * 128 - 64, z * 128 - 64), Voroce::Hash3DLowQuality).first;
		}
		const auto end = std::chrono::system_clock::now();
		const auto  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "3d opt: " << elapsed << std::endl;
	}

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto z = y;
			const auto t = x;
			volatile auto v = Voroce::Evaluate4DRef(glm::vec4(x * 128 - 64, y * 128 - 64, z * 128 - 64, t * 128 - 64), Voroce::Hash4DLowQuality).first;
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "4d ref: " << elapsed << std::endl;
	}

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto z = y;
			const auto t = x;
			volatile auto v = Voroce::Evaluate4DOpt(glm::vec4(x * 128 - 64, y * 128 - 64, z * 128 - 64, t * 128 - 64), Voroce::Hash4DLowQuality).first;
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "4d opt: " << elapsed << std::endl;
	}

	{
		std::mt19937 mt(15);
		std::uniform_real_distribution<double> uni(0.0f, 1.0f);
		for (auto L = 0; L < 32; ++L)
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
					const auto jitter = 1.0f;
					const auto ref = Voroce::Evaluate2DRef(glm::vec2(x * 128 - 64, y * 128 - 64), Voroce::Hash2DLowQuality, jitter);
					const auto opt = Voroce::Evaluate2DOpt(glm::vec2(x * 128 - 64, y * 128 - 64), Voroce::Hash2DLowQuality, jitter);
					if (ref.second != opt.second)
					{
						std::cout << "NG: ref:" << ref.first << " opt:" << opt.first << " ref:" << std::setprecision(16) << ref.second << " opt:" << opt.second << " jitter:" << jitter << std::endl;
					}
				});
		}
	}

	{
		std::mt19937 mt(15);
		std::uniform_real_distribution<double> uni(0.0f, 1.0f);
		for (auto L = 0; L < 32; ++L)
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
					const auto jitter = 1.0f;
					const auto ref = Voroce::Evaluate3DRef(glm::vec3(x * 128 - 64, y * 128 - 64, z * 128 - 64), Voroce::Hash3DLowQuality, jitter);
					const auto opt = Voroce::Evaluate3DOpt(glm::vec3(x * 128 - 64, y * 128 - 64, z * 128 - 64), Voroce::Hash3DLowQuality, jitter);
					if (ref.second != opt.second)
					{
						std::cout << "NG: ref:" << ref.first << " opt:" << opt.first << " ref:" << std::setprecision(16) << ref.second << " opt:" << opt.second << " jitter:" << jitter << std::endl;
					}
				});
		}
	}

	{
		// bug or numerical error
		std::mt19937 mt(15);
		std::uniform_real_distribution<double> uni(0.0f, 1.0f);
		for (auto L = 0; L < 32; ++L)
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
					const auto jitter = 1.0f;
					const auto ref = Voroce::Evaluate4DRef(glm::vec4(x * 128 - 64, y * 128 - 64, z * 128 - 64, t * 128 - 64), Voroce::Hash4DLowQuality, jitter);
					const auto opt = Voroce::Evaluate4DOpt(glm::vec4(x * 128 - 64, y * 128 - 64, z * 128 - 64, t * 128 - 64), Voroce::Hash4DLowQuality, jitter);
					if (ref.second != opt.second)
					{
						std::cout << "NG: ref:" << ref.first << " opt:" << opt.first << " ref:" << std::setprecision(16) << ref.second << " opt:" << opt.second << " jitter:" << jitter << std::endl;
					}
				});
		}
	}
}
