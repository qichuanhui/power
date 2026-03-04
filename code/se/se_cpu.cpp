#include "se.h" // 包含潮流计算核心函数
#include <iostream>
#include <cstdlib>
#include <sys/time.h>

int main(int argc, char *argv[]) {
    // 1. 参数检查
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <casename> <precision>" << std::endl;
        std::cerr << "Example: " << argv[0] << " Case10790 0.001" << std::endl;
        return EXIT_FAILURE;
    }

    // 2. 获取命令行参数
    std::string casename = argv[1];
    double precision = std::stod(argv[2]);
    
    std::cout << "Starting Power Flow Calculation..." << std::endl;
    std::cout << "Case: " << casename << ", Precision: " << precision << std::endl;
    
    // 3. 记录总开始时间
    struct timeval total_start, total_end;
    gettimeofday(&total_start, nullptr);
    
    // 4. 执行状态估计（选择单线程或多线程版本）
            // 单线程版本
    se_iter_time();  
    // 5. 记录总结束时间并计算总耗时
    gettimeofday(&total_end, nullptr);
    double total_time = (total_end.tv_sec - total_start.tv_sec) +
                       (total_end.tv_usec - total_start.tv_usec) / 1000000.0;
    std::cout << "data is:" << casename<<std::endl;
    std::cout << "\nPower Flow Calculation Completed!" << std::endl;
    std::cout << "Total Execution Time: " << total_time*100 << " seconds" << std::endl;
    
    return EXIT_SUCCESS;
}