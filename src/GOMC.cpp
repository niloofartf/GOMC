#include "InputFileReader.h"
#include <string>

std::string EnsembleSearch(std::string filename);

int main(int argc, char* argv[]) {

  std::string filename(argv[1]);
	std::string ensembleType = EnsembleSearch(filename);

	std::string Executable_To_Run = "./GOMC_CPU_";
	Executable_To_Run += ensembleType;
  Executable_To_Run += " ";
  Executable_To_Run += filename;

#ifdef _WIN32
	// TODO
#endif

#if defined(__linux__) || defined(__APPLE__)
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
