//
// Created by ubundrian on 21.06.16.
//

#include "Logger.h"

void Logger::openFile(std::string filename)
{
    std::string fileUrl = RESOURCES_PATH;
    fileUrl += "/logs/";
    fileUrl += filename;
    outputFile.open(fileUrl, std::ios::out);
}

void Logger::writeToFile(std::string text)
{
    if (outputFile.is_open()) {
        print("Write '" + text + "' to file");
        outputFile << text;
    } else {
        print("File has not been opened before writing", Mode::ERROR);
    }
}

void Logger::closeFile()
{
    outputFile.close();
}

void Logger::print(std::string text, Mode mode)
{
    if (mode & LOG)
    {
        std::cout << tabs << text << std::endl;
    } else if (mode & WARNING)
    {
        std::cout << tabs << "WARNING: " << text << std::endl;
    } else if (mode & ERROR)
    {
        std::cerr << tabs << text << std::endl;
    }
}

void Logger::tabIn()
{
    indent++;
    calcTab();
}

void Logger::tabOut()
{
    indent = std::max(0, indent-1);
    calcTab();
}

void Logger::changeTab(std::string newTab)
{
    tab = newTab;
}

void Logger::calcTab()
{
    tabs = "";
    for (int i = 0; i < indent; i++)
    {
        tabs += tab;
    }
}