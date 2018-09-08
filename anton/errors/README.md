Я проводил два сравнения: данные Антона из https://github.com/mryodo/nonlinSolver2sp/tree/master/matlab_main/data
с результатами своих расчётов на matlab (просто чтобы убедиться, что там лежит действительно то, что нужно) и
получившиеся на matlabе «свои» данные и данные вычислений на c++.

Если что, compare_among_matlab запускается примерно так:

./compare_among_matlab <my_N1> <my_N2> <anton_N1> <anton_N2> <low_x> <high_x> <low_y> <high_y> <ccto или hm> <dim(1, 2 или 3)>

А compare_cpp_vs_matlab так:

./compare_cpp_vs_matlab <ccto или hm> <dim(1, 2 или 3)> <папка с C++ данными> <папка с Matlab данными>


NB: ccto D1 — это вычисление main_ccto_d1_full.m (я убрал комментарии в оригинальном файле main_ccto_d1.m и сделал его похожим на остальные). Я не понял, где у Антона лежат результаты этого вычисления и не сравнивал их с тем, что на matlabe получилось у меня. (Если что, в «укороченной» версии ccto D1, где всего лишь один параметр, процентная разница между c++ и matlab'ом минимальная: см. график error_ccto_dim1_after_formatting.png.)
