/* This file is part of the C++ core of 'literanger'.
 *
 * literanger's C++ core was adapted from the C++ core of the 'ranger' package
 * for R Statistical Software <https://www.r-project.org>. The ranger C++ core
 * is Copyright (c) [2014-2018] [Marvin N. Wright] and distributed with MIT
 * license. literanger's C++ core is distributed with the same license, terms,
 * and permissions as ranger's C++ core.
 *
 * Copyright [2023] [Stephen Wade]
 *
 * This software may be modified and distributed under the terms of the MIT
 * license. You should have received a copy of the MIT license along with
 * literanger. If not, see <https://opensource.org/license/mit/>.
 */
#ifndef LITERANGER_UTILITY_INTERRUPT_H
#define LITERANGER_UTILITY_INTERRUPT_H

namespace literanger {

/** An interrupt check operator
 *
 * Defaults to never interrupted */
struct interruptor { virtual bool operator()() const; };

inline bool interruptor::operator()() const { return false; }


} /* namespace literanger */


#endif /* LITERANGER_UTILITY_INTERRUPT_H */

