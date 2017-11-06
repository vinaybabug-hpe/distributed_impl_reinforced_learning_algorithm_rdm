/*
 ============================================================================
 Name        : cclogger.h
 Author      : Vinay B Gavirangaswamy
 Created on	 : Nov 20, 2015
 Version     : 1.0
 Copyright   :  This program is free software: you can redistribute it and/or modify
    			it under the terms of the GNU General Public License as published by
    			the Free Software Foundation, either version 3 of the License, or
    			(at your option) any later version.

    			This program is distributed in the hope that it will be useful,
    			but WITHOUT ANY WARRANTY; without even the implied warranty of
    			MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    			GNU General Public License for more details.


    			You should have received a copy of the GNU General Public License
    			along with this program.  If not, see <http://www.gnu.org/licenses/>.
 Description : 
 ============================================================================
 */

#ifndef INC_COMMON_CCLOGGER_H_
#define INC_COMMON_CCLOGGER_H_

#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <common/wrapper.h>

class NullBuffer : public std::streambuf
{
public:
  int overflow(int c) { return c; }
};


enum LogSeverity { INFO, WARNING, ERROR, FATAL };

class Logger {
 public:
  Logger(LogSeverity ls, const std::string& file, int line,const std::string& method)
      : ls_(ls), file_(file), line_(line), method_(method)
  {


  }

  std::ostream& stream() const {

    return std::cerr << "["<<file_ << " | " << line_ <<" | " << method_<<"]: ";

  }
  ~Logger() {
    if (ls_ == FATAL) {

    }
  }
 private:
  LogSeverity ls_;
  std::string file_;
  int line_;
  std::string method_;



};

NullBuffer null_buffer;
std::ostream null_stream(&null_buffer);

#define LOG(ls)  PRINT_DEBUG==1? Logger(ls, __FILE__, __LINE__, __FUNCTION__).stream() : null_stream




#endif /* INC_COMMON_CCLOGGER_H_ */
