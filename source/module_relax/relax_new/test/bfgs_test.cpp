#include <gtest/gtest.h>
#include "../../mytest/bfgs.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/matrix3.h"
#include "module_base/vector3.h"

TEST(BFGSTest, InitializeTest) {
    bfgs optimizer;
    int size = 5;
    optimizer.initialize(size);

    EXPECT_EQ(optimizer.size, size);
    EXPECT_EQ(optimizer.alpha, 70);
    EXPECT_EQ(optimizer.maxstep, 100);
    for (int i = 0; i < 3 * size; ++i) {
        EXPECT_DOUBLE_EQ(optimizer.H[i][i], 70.0);
    }
}




// test PrepareStep 
TEST(BFGSTest, PrepareStepTest) {
    bfgs optimizer;
    int size = 3;
    optimizer.initialize(size);

    std::vector<std::vector<double>> test_pos(size, std::vector<double>(3, 0.0));
    std::vector<std::vector<double>> test_force(size, std::vector<double>(3, 0.5));

    optimizer.pos = test_pos;
    optimizer.force = test_force;
    optimizer.PrepareStep();
    for (const auto& row : optimizer.dpos) {
        for (const auto& elem : row) {
            EXPECT_NE(elem, 0.0);
        }
    }
}

// test IsRestrain 
TEST(BFGSTest, IsRestrainTest) {
    bfgs optimizer;
    int size = 3;
    optimizer.initialize(size);

    optimizer.dpos = {{1e-8, 1e-8, 1e-8}, {1e-8, 1e-8, 1e-8}, {1e-8, 1e-8, 1e-8}};
    EXPECT_TRUE(optimizer.IsRestrain());

    optimizer.dpos = {{0.001, 0.001, 0.001}, {0.001, 0.001, 0.001}, {0.001, 0.001, 0.001}};
    EXPECT_FALSE(optimizer.IsRestrain());
}

// test ReshapeMToV 
TEST(BFGSTest, ReshapeMToVTest) {
    bfgs optimizer;
    std::vector<std::vector<double>> matrix = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}}; 
    std::vector<double> expected = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}; 
    EXPECT_EQ(optimizer.ReshapeMToV(matrix), expected); 
    }
// test DotInVAndV  
TEST(BFGSTest, DotInVAndVTest) 
{ 
    bfgs optimizer; 
    std::vector<double> vec1 = {1.0, 2.0, 3.0}; 
    std::vector<double> vec2 = {4.0, 5.0, 6.0}; 
    double result = optimizer.DotInVAndV(vec1, vec2); 
    EXPECT_DOUBLE_EQ(result, 32.0); 
} 
