// APACHE LICENSE, VERSION 2.0
// copyright (c) shinji ogaki

#include <chrono>
#include <iostream>
#include <iomanip>

#include "../include/voroce.h"

int main()
{
	const auto golden_ratio = 1.61803398875;
	const auto N = 20000000;

	{
		const auto start = std::chrono::system_clock::now();
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			volatile auto v = Voroce::Evaluate2DRef(glm::vec2(x * 32 - 16, y * 32 - 16), Voroce::Hash2DLowQuality).first;
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
			volatile auto v = Voroce::Evaluate2DOpt(glm::vec2(x * 32 - 16, y * 32 - 16), Voroce::Hash2DLowQuality).first;
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
			const auto ref = Voroce::Evaluate2DRef(glm::vec2(x * 32 - 16, y * 32 - 16), Voroce::Hash2DLowQuality).first;
			const auto opt = Voroce::Evaluate2DOpt(glm::vec2(x * 32 - 16, y * 32 - 16), Voroce::Hash2DLowQuality).first;
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
			volatile auto v = Voroce::Evaluate3DRef(glm::vec3(x * 32 - 16, y * 32 - 16, z * 32 - 16), Voroce::Hash3DLowQuality).first;
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
			volatile auto v = Voroce::Evaluate3DOpt(glm::vec3(x * 32 - 16, y * 32 - 16, z * 32 - 16), Voroce::Hash3DLowQuality).first;
		}
		const auto end = std::chrono::system_clock::now();
		const auto  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		std::cout << "3d opt: " << elapsed << std::endl;
	}

	{
		// bug or numerical error
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto ref = Voroce::Evaluate3DRef(glm::vec3(x * 32 - 16, y * 32 - 16, y * 32 - 16), Voroce::Hash3DLowQuality).first;
			const auto opt = Voroce::Evaluate3DOpt(glm::vec3(x * 32 - 16, y * 32 - 16, y * 32 - 16), Voroce::Hash3DLowQuality).first;
			if (ref != opt)
			{
				std::cout << "NG: (" << x << "," << y << ") ref:" << ref << " opt:" << opt << std::endl;
			}
		}
	}

	{
		// bug or numerical error
		for (auto i = 0; i < N; ++i)
		{
			const auto x = (i + 0.5) / N;
			const auto y = (i * golden_ratio) - std::floor(i * golden_ratio);
			const auto ref = Voroce::Evaluate3DRef(glm::vec3(x * 32 - 16, y * 32 - 16, y * 32 - 16), Voroce::Hash3DLowQuality, 0.0f).first;
			const auto opt = Voroce::Evaluate3DOpt(glm::vec3(x * 32 - 16, y * 32 - 16, y * 32 - 16), Voroce::Hash3DLowQuality, 0.0f).first;
			if (ref != opt)
			{
				std::cout << "NG: (" << x << "," << y << ") ref:" << ref << " opt:" << opt << std::endl;
			}
		}
	}
}
