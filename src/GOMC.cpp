#include "InputFileReader.h"
#include <string>
#include <stdlib.h>
#include <iostream>

std::string EnsembleSearch(std::string filename);

int main(int argc, char* argv[])
{
  if(argc < 2) {
    std::cerr << "Error: Input file name is required!" << std::endl;
    std::cerr << "Example: GOMC in.conf" << std::endl;
    exit(0);
  }
  std::string filename(argv[1]);
	std::string ensembleType = EnsembleSearch(filename);

  if(ensembleType == "NAN") {
    std::cerr << "Error: Ensemble type is required in your input file!" << std::endl;
    std::cerr << "Example: Ensemble GCMC" << std::endl;
    exit(0);
  }

#ifdef _WIN32
	// TODO
#endif

#if defined(__linux__) || defined(__APPLE__)
	std::string Executable_To_Run = "./GOMC_CPU_";
	Executable_To_Run += ensembleType;
  Executable_To_Run += " ";
  Executable_To_Run += filename;
	system(Executable_To_Run.c_str());
#endif
	
	return 0;
}

std::string EnsembleSearch(std::string filename)
{
  std::vector<std::string> line;
  InputFileReader reader;
  reader.Open(filename);

	while (reader.readNextLine(line))
	{
		if (line.size() == 0)
			continue;
		if (line[0] == "Ensemble" || line[0] == "ENSEMBLE" || line[0] == "ensemble")
		{
			return line[1]; //return the word after ensemble
		}
	}
  return "NAN";
}
