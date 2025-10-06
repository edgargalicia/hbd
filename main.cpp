#include "logic.h"
#include <exception>

int main(int argc, char *argv[]) {
  try {
    Logic logic;
    logic.Init();
    // logic.Run();
    logic.Finish();
    return 0;
  }
  catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << '\n';
    return 1;
  }
}
