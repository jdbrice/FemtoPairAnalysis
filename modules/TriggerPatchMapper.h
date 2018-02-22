#ifndef TRIGGER_PATCH_MAPPER_H
#define TRIGGER_PATCH_MAPPER_H

#include "IObject.h"
using namespace jdb;


#include <map>
#include <utility>
#include <memory>
using namespace std;

class TriggerPatchMapper :public IObject
{
protected:
	// uses bl+1 (starts at 1)
	// and module +1 (starts at 1)
	// to make it easier to copy from pdf
	int patchMap[31][6];
public:
	virtual const char* classname() const { return "TriggerPatchMapper"; }
	TriggerPatchMapper() {

		for ( int iBL = 0; iBL < 31; iBL++ ){
			for ( int iM = 0; iM < 6; iM++ ){
				patchMap[iBL][iM] = -1;
			}
		}

		//////////start//////////
		patchMap[1 ][1] = 7; 
		patchMap[2 ][1] = 7; 
		patchMap[30][1] = 7; 
		patchMap[29][1] = 7; 
		patchMap[28][1] = 7; 

		patchMap[1 ][2] = 15;
		patchMap[2 ][2] = 15;
		patchMap[30][2] = 15;
		patchMap[29][2] = 15;
		patchMap[28][2] = 15;

		patchMap[1 ][3] = 6;
		patchMap[2 ][3] = 6;
		patchMap[30][3] = 6;
		patchMap[29][3] = 6;
		patchMap[28][3] = 6;

		patchMap[1 ][4] = 16;
		patchMap[2 ][4] = 16;
		patchMap[30][4] = 16;
		patchMap[29][4] = 16;
		patchMap[28][4] = 16;

		patchMap[1 ][5] = 8;
		patchMap[2 ][5] = 8;
		patchMap[30][5] = 8;
		patchMap[29][5] = 8;
		patchMap[28][5] = 8;
		////////stop////////////

		//////////start//////////
		patchMap[27][1] = 1; 
		patchMap[26][1] = 1; 
		patchMap[25][1] = 1; 
		patchMap[24][1] = 1; 

		patchMap[27][2] = 3;
		patchMap[26][2] = 3;
		patchMap[25][2] = 3;
		patchMap[24][2] = 3;

		patchMap[27][3] = 5;
		patchMap[26][3] = 5;
		patchMap[25][3] = 5;
		patchMap[24][3] = 5;

		patchMap[27][4] = 4;
		patchMap[26][4] = 4;
		patchMap[25][4] = 4;
		patchMap[24][4] = 4;

		patchMap[27][5] = 2;
		patchMap[26][5] = 2;
		patchMap[25][5] = 2;
		patchMap[24][5] = 2;
		////////stop////////////

		//////////start//////////
		patchMap[3][1] = 9; 
		patchMap[4][1] = 9; 
		patchMap[5][1] = 9; 
		patchMap[6][1] = 9; 
		patchMap[7][1] = 9; 

		patchMap[3][2] = 11;
		patchMap[4][2] = 11;
		patchMap[5][2] = 11;
		patchMap[6][2] = 11;
		patchMap[7][2] = 11;

		patchMap[3][3] = 13;
		patchMap[4][3] = 13;
		patchMap[5][3] = 13;
		patchMap[6][3] = 13;
		patchMap[7][3] = 13;

		patchMap[3][4] = 12;
		patchMap[4][4] = 12;
		patchMap[5][4] = 12;
		patchMap[6][4] = 12;
		patchMap[7][4] = 12;

		patchMap[3][5] = 10;
		patchMap[4][5] = 10;
		patchMap[5][5] = 10;
		patchMap[6][5] = 10;
		patchMap[7][5] = 10;
		////////stop////////////

		//////////start//////////
		patchMap[8 ][1] = 17; 
		patchMap[10][1] = 17; 
		patchMap[11][1] = 17; 
		patchMap[12][1] = 17; 

		patchMap[8 ][2] = 19;
		patchMap[10][2] = 19;
		patchMap[11][2] = 19;
		patchMap[12][2] = 19;

		patchMap[8 ][3] = 21;
		patchMap[10][3] = 21;
		patchMap[11][3] = 21;
		patchMap[12][3] = 21;

		patchMap[8 ][4] = 20;
		patchMap[10][4] = 20;
		patchMap[11][4] = 20;
		patchMap[12][4] = 20;

		patchMap[8 ][5] = 18;
		patchMap[10][5] = 18;
		patchMap[11][5] = 18;
		patchMap[12][5] = 18;
		////////stop////////////

		//////////start//////////
		patchMap[13][2] = 31;
		patchMap[14][2] = 31;
		patchMap[15][2] = 31;
		patchMap[16][2] = 31;
		patchMap[17][2] = 31;

		patchMap[13][3] = 22;
		patchMap[14][3] = 22;
		patchMap[15][3] = 22;
		patchMap[16][3] = 22;
		patchMap[17][3] = 22;

		patchMap[13][4] = 32;
		patchMap[14][4] = 32;
		patchMap[15][4] = 32;
		patchMap[16][4] = 32;
		patchMap[17][4] = 32;
		////////stop////////////


		//////////start//////////
		patchMap[22][1] = 25; 
		patchMap[21][1] = 25; 
		patchMap[20][1] = 25; 
		patchMap[19][1] = 25; 
		patchMap[18][1] = 25; 

		patchMap[22][2] = 27;
		patchMap[21][2] = 27;
		patchMap[20][2] = 27;
		patchMap[19][2] = 27;
		patchMap[18][2] = 27;

		patchMap[22][3] = 29;
		patchMap[21][3] = 29;
		patchMap[20][3] = 29;
		patchMap[19][3] = 29;
		patchMap[18][3] = 29;

		patchMap[22][4] = 28;
		patchMap[21][4] = 28;
		patchMap[20][4] = 28;
		patchMap[19][4] = 28;
		patchMap[18][4] = 28;

		patchMap[22][5] = 26;
		patchMap[21][5] = 26;
		patchMap[20][5] = 26;
		patchMap[19][5] = 26;
		patchMap[18][5] = 26;
		////////stop////////////

	}

	static shared_ptr<TriggerPatchMapper> instance;
	~TriggerPatchMapper() {}

	// 0 to 11
	static int cell ( int _channel ){
		return ( _channel % 60 ) % 12;
	}

	// 0 to 4
	static int module( int _channel ){
		return (_channel % 60 ) / 12;
	}
	// 0 to 29
	static int backleg( int _channel ){
		return _channel / 60;
	}

	static int etaStrip( int _channel ){
		int mod = TriggerPatchMapper::module( _channel ) + 5 * TriggerPatchMapper::backleg( _channel );
		return mod % 5 - 2;
	}

	int triggerPatch(  int _channel ){
		int bl = TriggerPatchMapper::backleg( _channel );
		int mod = TriggerPatchMapper::module( _channel );

		return patchMap[bl+1][mod+1];
	}

	static int findTriggerPatch( int _channel );
	
	int triggerPatch(  int _bl, int _mod ){
		return patchMap[_bl][_mod];
	}

	static int findTriggerPatch( int _bl, int _mod );

	
};




#endif