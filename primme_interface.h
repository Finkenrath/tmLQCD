/***********************************************************************
 *
 * Copyright (C) 2018 Jacob Finkenrath
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Interface for primme 
 *
 *******************************************************************************/

#ifndef PRIMME_INTERFACE_H_
#define PRIMME_INTERFACE_H_
#include "global.h"
#include "su3.h"
#include "solver/matrix_mult_typedef.h"
#include "solver/matrix_mult_typedef_nd.h"
#include "primme.h"

extern primme_params primme_tm;
void primme_tm_init(int op, int nev, int val, double prec);
void primme_tm_set_meth(int meth);
void primme_tm_reset(void);
void primme_tm_finalize(void);
int primme_tm_ev(void);
void primme_tmsvds_init(int aug);
void primme_tmsvds_set_meth(int id);
void primme_tmsvds_finalize(void);
int primme_tm_svds(void);


#endif /* PRIMME_INTERFACE_H_ */
