#include "rand.h"
#include <random>

template<typename T>
static constexpr T sqr(T a) {
    return a * a;
}

template<typename T>
static constexpr T power(T a, size_t n) {
    return n == 0 ? 1 : sqr(power(a, n / 2)) * (n % 2 == 0 ?  1 : a);
}

int random(int begin, int end)
{
    static std::mt19937 gen { std::random_device{}() };
//    static std::mt19937 gen{0};
    static std::uniform_int_distribution<int> dist(0, std::numeric_limits<int>::max());
    return (dist(gen) % (end - begin + 1)) + begin;
}

double random(double begin, double end, unsigned int precision)
{
    const int divisor = power(10, (precision > 20) ? 20 : precision);
    int i_begin = static_cast<int>(begin * divisor);
    int i_end = static_cast<int>(end * divisor);
    return random(i_begin, i_end) / static_cast<double>(divisor);
}
