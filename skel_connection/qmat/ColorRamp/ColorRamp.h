//********************************************
// ColorRamp.h
//********************************************
// pierre.alliez@cnet.francetelecom.fr
// Created : 19/05/98
// Modified : 19/05/98
//********************************************

#ifndef _COLOR_RAMP_
#define _COLOR_RAMP_


// Datas : 
// Red Green Blue IsNode (0/1)

class CColorRamp
{
private :
	int m_Color[4][256];
	int m_Node[256];
	int  m_Size;
	int  m_NbNode;

public :

	// Constructor
	CColorRamp();
	~CColorRamp();

	// Datas
	int GetSize() { return m_Size; }
private:
	int Red(int index) { if(index < 0) index = 0; if(index > 255) index = 255; return m_Color[0][index]; }
	int Green(int index) { if(index < 0) index = 0; if(index > 255) index = 255; return m_Color[1][index]; }
	int Blue(int index) { if(index < 0) index = 0; if(index > 255) index = 255; return m_Color[2][index]; }
public:
	void RedGreenBlue(int index, float mc[3]);
	

	// Misc
	int Build();
	void BuildDefault();
	int BuildNodes();
	void ResetNodes();

	void BuildRainbow();

};

#endif // _COLOR_RAMP_
