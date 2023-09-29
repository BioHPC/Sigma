#ifndef PARSE_CONFIG_H
#define PARSE_CONFIG_H

//=============================================================================
// parse_config.h
//   : parse_config.h is the header file for parssing config file
//
// Author: Tae-Hyuk (Ted) Ahn
// Created: 04/01/2013
// Modified:
//
// Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). All rights reserved.
//=============================================================================


#include "utils.h"


//=============================================================================
// Class of parse config
//=============================================================================
class ParseConfig
{
public:
    // Default constructor
    ParseConfig(const std::string &config_path);

    // Default destructor
    ~ParseConfig();

    // Extarct keys and values, then save to map
    void extractConfigFile(const std::string &config_path);

    // Check valid line
    bool validLine(const std::string &line);

    // Insert key->value to config_map
    void insertKeyValueMap(const std::string &line, const std::string &section_name);

    // Extract key from line
    void extractKey(std::string &key, size_t const &sepPos, const std::string &line);

    // Extract value from line
    void extractValue(std::string &value, size_t const &sepPos, const std::string &line);

    // Check key exists in config_map
    bool keyExists(const std::string &key);

    // Get value of key
    std::string getValueOfKey(const std::string &key);

private:
    // Private default constructor
    ParseConfig() { }

    // Config config_path
    static std::string config_path;

    // Config map
    std::map<std::string, std::string> *config_map;

};

#endif
