bool passed_all = true;

#define IS_TRUE(x) {if (!(x)) {std::cout << __FUNCTION__ << " failed on line " << __LINE__ << std::endl; passed_all = false;}}
#define IS_FALSE(x) {IS_TRUE(!x);}