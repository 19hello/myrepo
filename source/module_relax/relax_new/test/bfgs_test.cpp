/*#include <gtest/gtest.h>
#include "../bfgs.h"
#include "module_hamilt_pw/hamilt_pwdft/global.h"
#include "module_base/matrix3.h"
#include "module_base/vector3.h"

TEST(BFGSTest, InitializeTest) {
    BFGS optimizer;
    int size = 5;
    optimizer.init_relax(size);

    EXPECT_EQ(optimizer.size, size);
    EXPECT_EQ(optimizer.alpha, 70);
    EXPECT_EQ(optimizer.maxstep, 100);
    for (int i = 0; i < 3 * size; ++i) {
        EXPECT_DOUBLE_EQ(optimizer.H[i][i], 70.0);
    }
}




// test PrepareStep 
TEST(BFGSTest, PrepareStepTest) {
    BFGS optimizer;
    int size = 3;
    optimizer.init_relax(size);

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
    BFGS optimizer;
    int size = 3;
    optimizer.init_relax(size);

    optimizer.dpos = {{1e-8, 1e-8, 1e-8}, {1e-8, 1e-8, 1e-8}, {1e-8, 1e-8, 1e-8}};
    EXPECT_TRUE(optimizer.IsRestrain());

    optimizer.dpos = {{0.001, 0.001, 0.001}, {0.001, 0.001, 0.001}, {0.001, 0.001, 0.001}};
    EXPECT_FALSE(optimizer.IsRestrain());
}*/


