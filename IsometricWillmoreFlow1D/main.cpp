#include "IsometricWillmoreFlow1D.h"

int main(int argc, char *argv[])
{
	if (argc < 2) return -1;

	auto flow = new IsometricWillmoreFlow1D(new Mesh(argv[1]));

	while (true)
	{
		flow->integrate(0.01);
		flow->mesh->center();
	}

    return 0;
}
