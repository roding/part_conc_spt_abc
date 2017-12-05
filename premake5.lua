workspace "part_conc_accel"
	platforms { "native", "x64", "x32" }
	configurations { "debug", "release" }

	language "C++"
	cppdialect "C++11" -- target NVCC (8.0) doesn't support C++14 yet. :-(

	flags "NoPCH"
	flags "MultiProcessorCompile"

	vectorextensions "AVX"

	CUDAgencode = {
		"arch=compute_52,code=sm_52",
		--"arch=compute_60,code=sm_60",
		"arch=compute_61,code=sm_61"
	};
	CUDAflags = {
		"--expt-relaxed-constexpr",
		"-Xptxas='-v'"
	};

	objdir "build/%{cfg.buildcfg}-%{cfg.platform}-%{cfg.toolset}"
	targetsuffix "-%{cfg.buildcfg}-%{cfg.platform}-%{cfg.toolset}"

	newoption {
		trigger = "toolset",
		description = "Select toolset on Linux / MacOS",
		allowed = {
			{ "gcc", "Use GCC" },
			{ "clang", "Use Clang" }
		}
	};

	-- CUDA/nvcc builds
	cudaCodeStr = "";
	for _,gc in ipairs(CUDAgencode) do
		cudaCodeStr = cudaCodeStr .. "-gencode " .. gc .. " ";
	end
	cudaFlagStr = "";
	for _,flag in ipairs(CUDAflags) do
		cudaFlagStr = cudaFlagStr .. flag .. " ";
	end

	filter "files:**.cu"
		-- Not sure if this is the proper way to do this...
		local makeIPATH = [[%{(function()
			local ret = "";
			for _,path in ipairs(cfg.includedirs) do
				local istr = cfg.includedirs[path] or path;
				ret = ret .. " -I" .. istr;
			end
			return ret .. " ";
		end)()}]];
		local makeDEFS = [[%{(function()
			local ret = "";
			for _,def in ipairs(cfg.defines) do
				local dstr = cfg.defines[def] or def;
				ret = ret .. " -D" .. dstr;
			end
			return ret .. " ";
		end)()}]];
		local makeCCOM = [[%{(function()
			local ret = "-Xcompiler '";
			for _,opt in ipairs(cfg.buildoptions) do
				local ostr = cfg.buildoptions[opt] or opt;
				ret = ret .. ostr .. " ";
			end
			ret = ret .. "' ";
			return ret;
		end)()}]];
	
		buildmessage "CUDA: %{file.relpath}"
		buildcommands {
			"mkdir -p \"%{cfg.objdir}\"",
			"nvcc -std c++11 -g -O2 " .. cudaFlagStr .. cudaCodeStr .. makeIPATH .. makeDEFS .. makeCCOM .. " -c -o \"%{cfg.objdir}/%{file.basename}.o\" \"%{file.relpath}\""
		}
		buildoutputs { "%{cfg.objdir}/%{file.basename}.o" }
	filter "*"

	-- system libraries
	filter "system:linux"
		includedirs "/opt/cuda/include"
		includedirs "/usr/local/cuda-8.0/include"
		libdirs "/opt/cuda/lib64/"
		libdirs "/usr/local/cuda-8.0/lib64"
		links "cudart_static"
		links "dl"
		links "rt"
	filter "system:windows"
	
	filter "*"

	-- premake-workaround empty "toolset": always set it
	filter "system:linux or system:macos"
		toolset( _OPTIONS["toolset"] or "gcc" );
	filter "system:windows"
		toolset( "msc" );
	filter "*"

	-- default compiler flags
	filter "toolset:gcc or toolset:clang"
		linkoptions { "-pthread" }
		buildoptions { 
			"-march=native", 
			"-Wall", "-Wextra", 
			"-ffast-math", 
			"-pthread",
			"-Wno-missing-field-initializers"
		}

	filter "toolset:msc"
		defines { "_CRT_SECURE_NO_WARNINGS=1" }

	filter "action:vs2015"
		buildoptions { "/utf-8" }
	
	filter "action:vs2017"
		buildoptions { "/utf-8" }
		buildoptions { "/std:c++latest" }
	
	filter "*"

	-- default outputs
	filter "kind:StaticLib"
		targetdir "lib/"

	filter "kind:ConsoleApp"
		targetdir "bin/"
		targetextension ".exe"
	
	filter "*"

	--configurations
	configuration "debug"
		symbols "On" -- new, but doesn't work?
		optimize "On"
		defines { "_DEBUG=1" }

	configuration "release"
		optimize "On"
		defines { "NDEBUG=1" }

	configuration "*"


	-- projects
	project "support"
		kind "StaticLib"

		includedirs "external"

		files{ "src/gpu/support/**.cpp" }

	project "shared"
		kind "StaticLib"

		includedirs "external"

		files{ "src/gpu/shared/**.cpp" }

	project "cusim"
		kind "StaticLib"

		includedirs "external"

		files{ "src/gpu/cusim/*.cu" }

	project "accel"
		kind "ConsoleApp"

		includedirs "external"
		includedirs "src"

		files{ "src/gpu/*.cpp", "external/**.cpp", "src/gpu/*.cu" }

		links "cusim"
		links "shared"
		links "support"

	project "compare"
		kind "ConsoleApp"

		includedirs "external"
		includedirs "src"

		files{ "src/compare_output/**.cpp", "external/**.cpp" }
		links "shared"

--EOF
