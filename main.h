//#ifdef _WIN32
//#include <tchar.h>
//#endif
//
//  main.h
//  FEM_A02
//
//  Created by Wilson Wesley Wutzow and adapted by Paulo Henrique de Freitas Meirelles on 09/12/11.
//  Copyright 2011 Wutzow. All rights reserved.
//
#ifndef HmainH
#define HmainH
//--------------- Bibliotecas ----------------------------
#include <ostream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <math.h>
#include <string>
#include <cstring>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>

#include "CCO.h"
#include "ALG.h"
#include "NOS.h"
#include "ELM.h"
#include "MAT.h"
#include "QUA.h"

void Gera_VTK_69(CCo& Co, tvNo& No, tvEl& El, const int& ite);  // Gerador de saida VTK com endereçamento em relação as classes construtoras;


#endif
