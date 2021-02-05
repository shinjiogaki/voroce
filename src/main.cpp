// APACHE LICENSE, VERSION 2.0
// copyright (c) shinji ogaki

#include <chrono>
#include <iostream>
#include <iomanip>

#include "voroce.h"

using namespace voroce;

int main()
{
	const auto golden_ratio = 1.61803398875f;
	const auto N = 10498646;

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5f) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			volatile float v = Voroce::Evaluate2DRef(glm::vec2(x * 2 - 1, y * 2 - 1));
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "2d ref: " << elapsed << std::endl;
	}

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5f) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			volatile float v = Voroce::Evaluate2DOpt(glm::vec2(x * 2 - 1, y * 2 - 1));
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "2d opt: " << elapsed << std::endl;
	}

	{
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5f) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto ref = Voroce::Evaluate2DRef(glm::vec2(x * 2 - 1, y * 2 - 1));
			const auto opt = Voroce::Evaluate2DOpt(glm::vec2(x * 2 - 1, y * 2 - 1));
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
			const auto x = (i + 0.5f) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto z = y;
			volatile float v = Voroce::Evaluate3DRef(glm::vec3(x * 2 - 1, y * 2 - 1, z * 2 - 1));
		}
		const auto end = std::chrono::system_clock::now();
		const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "3d ref: " << elapsed << std::endl;
	}

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5f) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto z = y;
			volatile float v = Voroce::Evaluate3DOpt(glm::vec3(x * 2 - 1, y * 2 - 1, z * 2 - 1));
		}
		const auto end = std::chrono::system_clock::now();
		const auto  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "3d opt: " << elapsed << std::endl;
	}

	{
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5f) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto ref = Voroce::Evaluate3DRef(glm::vec3(x * 2 - 1, y * 2 - 1, y * 2 - 1));
			const auto opt = Voroce::Evaluate3DOpt(glm::vec3(x * 2 - 1, y * 2 - 1, y * 2 - 1));
			if (ref != opt)
			{
				std::cout << "NG: (" << x << "," << y << ") ref:" << ref << " opt:" << opt << std::endl;
			}
		}
	}
}
