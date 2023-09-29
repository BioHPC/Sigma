//=============================================================================
// parse_config.cpp
//   : parse_config.cpp is cpp file for parsing config file
//
// Author: Tae-Hyuk (Ted) Ahn
// Created: 04/01/2013
// Modified:
//
// Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). All rights reserved.
//=============================================================================


#include "utils.h"
#include "parse_config.h"


//=============================================================================
// Default constructor
//=============================================================================
ParseConfig::ParseConfig(const std::string &config_path)
{   
    config_map = new std::map<std::string, std::string>;
    extractConfigFile(config_path);
}


//=============================================================================
// Default destructor
//=============================================================================
ParseConfig::~ParseConfig()
{   
    delete config_map;
}


//=============================================================================
// Extarct keys and values, then save to map 
//=============================================================================
void ParseConfig::extractConfigFile(const std::string &config_path)
{
    // open file
    std::ifstream file;
    file.open(config_path.c_str());

    // string 
    std::string line, section_name;

    // line number
    unsigned int line_number = 0;

    // load config file
    while (std::getline(file, line))
    {

        // trim invisible spaces
        line.erase(line.find_last_not_of(" \n\r\t")+1);

        // line number
        line_number++;

        // empty line skip
        if (line.empty()) continue;
        // white space skip
        if (line.find_first_not_of(' ') == line.npos) continue;

        // comment line skip
        // if (line.find('#') != line.npos) line.erase(line.find('#'));
        if (line.at(0) == '#') continue;

        // for section name
        if ((line.at(0) == '[') && (line.at(line.length() -1) == ']')) {
            section_name = line;
        }
        // consider config key and value    
        else {
            // check if line has =
            if (line.find('=') == line.npos) {
                Utils::exitWithError("*** Error: Couldn't find separator in the config file on line: " + Utils::unsignedIntToString(line_number));
            }

            // check line has valid key and value
            if (!validLine(line)) {
                Utils::exitWithError("*** Error: Bad format in the config file on line: " + Utils::unsignedIntToString(line_number));
            }

            // insert key value map
            insertKeyValueMap(line, section_name);
        }
    }
        
    file.close();
}


//=============================================================================
// Check valid line 
//=============================================================================
bool ParseConfig::validLine(const std::string &line)
{
    std::string temp = line;
    temp.erase(0, temp.find_first_not_of("\t "));
    if (temp[0] == '=')
        return false;

    for (size_t i = temp.find('=') + 1; i < temp.length(); i++)
        if (temp[i] != ' ')
            return true;

    return false;
}


//=============================================================================
// Insert key->value to config_map 
//=============================================================================
void ParseConfig::insertKeyValueMap(const std::string &line, const std::string &section_name)
{
    std::string temp = line;
    temp.erase(0, temp.find_first_not_of("\t "));
    size_t sepPos = temp.find('=');

    std::string key, value;
    extractKey(key, sepPos, temp);
    key = section_name + key;
    extractValue(value, sepPos, temp);

    if (!keyExists(key)) {
        (*config_map).insert(std::pair<std::string, std::string>(key, value));
    }
    else {
        Utils::exitWithError("*** 1. Error: Can only have unique key names!");
    }
}


//=============================================================================
// Extract key from line 
//=============================================================================
void ParseConfig::extractKey(std::string &key, size_t const &sepPos, const std::string &line)
{
    key = line.substr(0, sepPos);
    if (key.find('\t') != line.npos || key.find(' ') != line.npos)
        key.erase(key.find_first_of("\t "));
}


//=============================================================================
// Extract value from line 
//=============================================================================
void ParseConfig::extractValue(std::string &value, size_t const &sepPos, const std::string &line)
{
    value = line.substr(sepPos + 1);
    value.erase(0, value.find_first_not_of("\t "));
    value.erase(value.find_last_not_of("\t ") + 1);
}


//=============================================================================
// Check key exists in config_map 
//=============================================================================
bool ParseConfig::keyExists(const std::string &key)
{
    return (*config_map).find(key) != (*config_map).end();
}


//=============================================================================
// Get value of key in config_map 
//=============================================================================
std::string ParseConfig::getValueOfKey(const std::string &key)
{
    if (!keyExists(key))
        Utils::exitWithError("*** Error: key " + key + " does not exist!");
    return (*config_map).find(key)->second;
}

