#include "aux_func.h"
#include <string>
using namespace std;

int main (int argc, char *argv[]) {

	genConsistDenseSystems();

	genConsistDenseSystemsNorm();

	genCoherentDenseSystems();

	return 0;
}