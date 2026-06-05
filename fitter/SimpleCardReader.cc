#include "SimpleCardReader.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TSystem.h"
#include <iostream>

const unsigned int SimpleCardReader::kMaxIncludeDepth = 5;

SimpleCardReader::SimpleCardReader() {
    // default constructor. do not call ReadCard
   kVerbosity = false;
   commentChars.push_back('#');
}

SimpleCardReader::SimpleCardReader( const char * cname , bool k)
{
   kVerbosity = k ;
   std::string CardName( cname );
   commentChars.push_back('#');
   ReadCard( CardName );
}

SimpleCardReader::SimpleCardReader( std::string CardName , bool k)
{
   kVerbosity = k ;
   commentChars.push_back('#');
   ReadCard( CardName );
}

void SimpleCardReader::ReadCard( std::string CardName, unsigned int depth )
{
    if (gSystem->AccessPathName(CardName.c_str(), kFileExists)) {
        fprintf(stderr, "[Error] SimpleCardReader::ReadCard %s does not exist.\n", CardName.c_str());
        abort();
    }
    if (depth > kMaxIncludeDepth) {
        fprintf(stderr, "[Error] SimpleCardReader::ReadCard maximum include depth exceeded. Do you have a loop?\n");
        abort();
    }

   std::ifstream InFile;
   InFile.open( CardName.c_str(), std::ifstream::in );
	
   if( !InFile.is_open()  )
   {
      std::cout << "SimpleCardReader::ReadCard " << CardName << " is not valid. " << std::endl; 
   }
   std::cout << "SimpleCardReader::ReadCard reading input from " << CardName << std::endl;

   std::string line;
   std::string prevline;
   std::string Value;
   std::string Name;
   std::string::size_type loc;

   int i = 1;
   int LineCount = 0;
   int ntokens_space;
   int ntokens_comma;
   int ntokens_string;
   bool kComment;
   std::vector < std::string > tokens;      
   std::vector < std::string > tokensTmp;      


   // Reset the list of keys
   Keys.clear();
   while( !InFile.eof() )
   {  
      getline( InFile, line);

      /// search for keys that span multiple lines
      std::stringstream ss;
      loc = line.find( '\\' , 0 ) ;
      ss << line ;
      prevline = line;
      while ( loc == prevline.length()-1 && loc != std::string::npos )
      {
         //std::cout <<" Processing multi-line at: " << line << "loc: " << loc <<  std::endl;
         getline( InFile, line);

         // allow comments mid-block, but skip them
         if ( CheckIfComment( line ) )  continue ;

         ss << line ;
         loc = line.find( '\\' ,0 ); 
         prevline = line;
      }  
      line = ss.str();
      std::string bs = "\\";
      Erase( line , bs.c_str()  );


      LineCount ++; 
      if ( LineCount > 10000 ) 
      { 
         std::cerr << "SimpleCardReader::ReadCard Too many lines in card file, 10000" << std::endl;
         abort();
      }

      // Skip lines with leading comment chars ('#' by default)
      if ( CheckIfComment( line ) ) continue;

      // turn any tabs into spaces
      Replace( line, "\t", " " );

      // search for simple variable definitions 
      ntokens_comma  = Tokenize( line , "," , tokensTmp ); 
      ntokens_space  = Tokenize( line , " " , tokensTmp ); 

 
      tokens.clear();
      Name  = tokensTmp[0];
      Value = "";
      for( unsigned k = 1 ; k < ntokens_space ; k++ )
         if( tokensTmp[k].length() > 0 &&  tokensTmp[k] != " " )
            tokens.push_back( tokensTmp[k] );
       
      // second loop so we do not place a spurious comma at the 
      // of sequences with a tail of spaces      
      for( unsigned k = 0 ; k < tokens.size() ; k++ )
      {
         Value += tokens[k];
         if( ntokens_comma == 1 )
           if( k < tokens.size()-1)
             Value += ","; 
      }

      // remove extra spaces
      Erase( Name , " " );
      Erase( Value, " " );

      // we are just strings now so forget about quotes 
      Erase( Value, "\"" ); 
   
      ExpandEnvironmentVariables( Value );

      // remove braces for array variables
      Erase( Value, "{" ); 
      Erase( Value, "}" );


      // import keys and values from another card
      if ( Name == "Include" ) {
        printf("SimpleCardReader::ReadCard Including keys from %s\n", Value.c_str());
        SimpleCardReader c;
        c.ReadCard(Value, depth + 1);
        std::vector<std::string> inc_keys;
        c.GetListOfKeys(inc_keys);
        for (const auto&  k : inc_keys) {
            // return string of key, regardless of what it actually is
            std::string v;
            c.GetKey(k, v);
            StringKeys[k] = v;
            if( tokens.size() == 1 )
              Keys[ k ] = atof( v.c_str() );
        }
      }

      StringKeys[ Name ] = Value;
    
      // here be kludge, fix it .
      if( tokens.size() == 1 )
         Keys[ Name ] = atof( Value.c_str() );
         

   }// end of loop over lines in card file

   if( kVerbosity )
      for( iSMap = StringKeys.begin() ; iSMap != StringKeys.end() ; iSMap++ )
          std::cout << iSMap->first << " [" << iSMap->second <<"]"<< std::endl;

}

void SimpleCardReader::KeyFindFailure( std::string key , const char * type , void * ptr )
{
  std::cout << "SimpleCardReader::GetKey \"" << key << "\" for type \"" << type << "\" was not found."
            << " variable at address " << ptr << " is unchanged." << std::endl; 
}


bool SimpleCardReader::FindKey( std::string SearchKey )
{
   return ( Keys.count( SearchKey ) > 0 ? true : false );
}


bool SimpleCardReader::GetKey( std::string SearchKey, double & V )
{
   if ( FindKey(SearchKey) )
   {
       V = Keys[ SearchKey ];
       return  true;
   }

   if( kVerbosity ) KeyFindFailure( SearchKey , "double" ,  &V );
   return false;
}


bool SimpleCardReader::GetKey( std::string SearchKey, float  & V )
{
   if ( FindKey(SearchKey) )
   {
       V = (float) Keys[ SearchKey ];
       return  true;
   }

   if( kVerbosity ) KeyFindFailure( SearchKey , "float" ,  &V );
   return false;
}
bool SimpleCardReader::GetKey( std::string SearchKey, int    & V )
{
   if ( FindKey(SearchKey) )
   {
       V = (int) Keys[ SearchKey ];
       return  true;
   }

   if( kVerbosity ) KeyFindFailure( SearchKey , "int" ,  &V );
   return false;
}
bool SimpleCardReader::GetKey( std::string SearchKey, bool   & V )
{
   if ( FindKey(SearchKey) )
   {
       V = Keys[ SearchKey ] > 0 ? true : false;
       return  true;
   }

   if( kVerbosity ) KeyFindFailure( SearchKey , "bool" ,  &V );
   return false;
}

bool SimpleCardReader::GetKey( std::string SearchKey, std::string & V )
{

   if ( StringKeys.count(SearchKey) > 0  )
   {
       V = StringKeys[ SearchKey ];
       return  true;
   }

   if( kVerbosity ) KeyFindFailure( SearchKey , "string" ,  &V );
   return false;
}
 
bool SimpleCardReader::GetKey( const char * s , std::string & V )  
{
   std::string SearchKey( s );
   return GetKey( SearchKey , V );
}

/// will automatically cut spaces and {} out std::string
int SimpleCardReader::ParseArray( std::string Value, const char *  D, double ** Out )
{

   double dValue;
   double  * tmp;

   std::vector< std::string > tokens;
   
   int ntokens = Tokenize( Value , D , tokens );   

   // add two spaces so that [0] stores number of bins
   ntokens+= 2;

   tmp = new double [ ntokens ];
   tmp[0] = (double)  ntokens  ;

   for( unsigned int i = 0 ; i < tokens.size(); i ++ )
     tmp[i+1] = atof( tokens[i].c_str() ); 
     

   *Out = &tmp[0];

   // the number of elements found
   return ntokens-2;

}


void SimpleCardReader::GetListOfKeys( std::vector< std::string >& v )
{
   v.clear();
   for( iSMap = StringKeys.begin() ; iSMap != StringKeys.end() ; iSMap++ )
      v.push_back( iSMap->first );

   return;
}

void SimpleCardReader::Erase( std::string & source, const std::string& kill, std::string::size_type l0 )
{

   std::string::size_type loc0 = source.find( kill, l0 );
  
   while( loc0 != std::string::npos )
   {
      source.erase( loc0, 1);
      loc0 = source.find(kill, loc0);
   }


}



unsigned int SimpleCardReader::Tokenize(const std::string& source,
                                        const std::string& delimiters,
                                        std::vector<std::string>& tokens)
{
    std::string::size_type prev_loc0 = 0;
    std::string::size_type loc0 = 0;
    std::string::size_type loc1 = 0;
    std::string sub;
    unsigned int ntokens = 0;

    tokens.clear();

    loc0 = source.find_first_of(delimiters, loc0);
    while (loc0 != std::string::npos)
    {

        sub = source.substr(prev_loc0, loc0 - prev_loc0);

        tokens.push_back( sub ) ;
        ntokens++;


        loc0++;
        prev_loc0 = loc0;
        loc0 = source.find_first_of(delimiters, loc0);
    }


    if (prev_loc0 < source.length())
    {
        tokens.push_back(source.substr(prev_loc0));
        ntokens++;
    }

    return ntokens;
}

void SimpleCardReader::Replace( std::string & source, const std::string & kill , 
                                                const std::string & rep  , std::string::size_type l0 )
{

   std::string::size_type loc0 = source.find( kill, l0 );
  
   while( loc0 != std::string::npos )
   {
      source.replace( loc0, kill.length(), rep );
      loc0 = source.find(kill, loc0);
   }

}

#include <regex>
#include <cstdlib>
void SimpleCardReader::ExpandEnvironmentVariables (std::string & val) {
    // static const std::regex env_re{R"--(\$\{([^}]+)\})--"};
    static const std::regex env_re{"\\$\\{([^}]+)\\}"};
    std::smatch match;
    while (std::regex_search(val, match, env_re)) {
        auto const from = match[0];
        auto const var_name = match[1].str().c_str();
        val.replace(from.first, from.second, std::getenv(var_name));
    }
}

bool SimpleCardReader::CheckIfComment( std::string line )
{

   bool kComment = false;
   std::string::size_type lloc;
   for( unsigned k = 0 ; k < commentChars.size() ; k++ )
   {
     lloc = line.find( commentChars[k] , 0 );
     if( lloc == 0 || line.empty() )
     {
       if( kVerbosity )
         std::cout << " Omitting Line as Comment: " << line << " (empty: " << line.empty() << ")" <<  std::endl;
       kComment = true ;
     }
   }

  return kComment;

}