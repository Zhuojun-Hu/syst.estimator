#ifndef SimpleCardReader_H
#define SimpleCardReader_H

#include <cstdlib>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "TString.h"


class SimpleCardReader
{
   public:
     
      SimpleCardReader();
      SimpleCardReader( const char *  , bool k = false );
      SimpleCardReader( std::string , bool k = false );

      void ReadCard( std::string, unsigned int depth = 0 );

      void KeyFindFailure( std::string key , const char * type , void * ptr );
      bool GetKey( std::string, double & );
      bool GetKey( std::string, float  & );
      bool GetKey( std::string, int    & );
      bool GetKey( std::string, bool   & );
      bool GetKey( std::string, std::string & );
      bool GetKey( const char*, std::string & );
  
      bool FindKey( std::string );

      int ParseArray( std::string , const char * , double ** );

      void GetListOfKeys( std::vector< std::string >& v );

      unsigned int Tokenize(const std::string & source,
                            const std::string & delimiters,
                            std::vector<std::string> & tokens);

      bool CheckIfComment( std::string line );
   private:      
      void Replace( std::string & , const std::string & , 
                                    const std::string & ,
                                    std::string::size_type l0 = 0 );
      void Erase( std::string & , const std::string & , std::string::size_type l0 = 0 );
      void ExpandEnvironmentVariables (std::string & );


      std::map< std::string , double >  Keys;
      std::map< std::string , double >::iterator iMap;

      std::map< std::string , std::string >  StringKeys;
      std::map< std::string , std::string >::iterator iSMap;
      
      std::vector< char > commentChars;
    
      bool kVerbosity;

      // for including other card files. Prevent too many recursive includes
      static const unsigned int kMaxIncludeDepth;

};

#endif



