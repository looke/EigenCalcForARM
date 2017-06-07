/*
 * gtest_main.cpp
 *
 *  Created on: 2017Äê6ÔÂ7ÈÕ
 *      Author: looke
 */

#include <stdio.h>
#include "..\gtest_src\gtest\gtest.h"

GTEST_API_ int main(int argc, char **argv)
{
	printf("Running main() from gtest_main.cc\n");
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

