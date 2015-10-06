-- fluid.lua
workspace "Fluid"
   configurations { "Debug", "Release" }
   location "build/Fluid"
   architecture "x86_64"

project "Fluid"
   kind "ConsoleApp"
   language "C++"
   targetdir "../bin/%{cfg.buildcfg}"

   files { "../src/**.h", "../src/*.cpp" }

   includedirs { "../3rdparty/GLFW/include"
   , "../3rdparty/GLEW/include"
   , "../3rdparty/glm/include"
   , "../3rdparty/boost/include" }

   libdirs { "../3rdparty/GLFW/lib"
   , "../3rdparty/GLEW/lib"
   , "../3rdparty/boost/lib" }

   links { "glfw3", "glew32", "glu32", "opengl32" }

   filter "configurations:Debug"
      defines { "DEBUG" }
      flags { "Symbols" }

   filter "configurations:Release"
      defines { "NDEBUG" }
      optimize "On"