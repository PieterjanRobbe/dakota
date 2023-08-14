/*  _______________________________________________________________________

    DAKOTA: Design Analysis Kit for Optimization and Terascale Applications
    Copyright 2014-2023
    National Technology & Engineering Solutions of Sandia, LLC (NTESS).
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Dakota directory.
    _______________________________________________________________________ */

#ifndef DAKOTA_BIT_UTILS_H
#define DAKOTA_BIT_UTILS_H

#include "dakota_data_types.hpp"

namespace Dakota {

/// Bit operations are taken from the Stanford Bithack webpage
/// http://www.graphics.stanford.edu/~seander/bithacks.html

/// Count consecutive trailing zero bits
inline unsigned count_consecutive_trailing_zero_bits(UInt32 v)
{
  unsigned c = 32;
  v &= -signed(v);
  if (v) c--;
  if (v & 0x0000FFFF) c -= 16;
  if (v & 0x00FF00FF) c -= 8;
  if (v & 0x0F0F0F0F) c -= 4;
  if (v & 0x33333333) c -= 2;
  if (v & 0x55555555) c -= 1;
  return c;
}

/// Reverse bits of unsigned 32 bit integer
inline UInt32 bitreverse(UInt32 k)
{
  UInt32 v = k;
  v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);  // swap odd and even bits
  v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);  // swap consecutive pairs
  v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);  // swap nibbles ...
  v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);  // swap bytes
  v = ( v >> 16             ) | ( v               << 16); // swap 2-byte long pairs
  return v;
}

/// Reverse bits of unsigned 64 bit integer
inline UInt64 bitreverse(UInt64 v)
{
  return (UInt64(bitreverse(UInt32(v))) << 32)
    | UInt64(bitreverse(UInt32(v >> 32)));
}

/// Convert binary to Gray code order
inline UInt64 binary2gray(UInt64 v)
{
  return v ^ (v >> 1);
}


} // namespace Dakota

#endif // DAKOTA_BIT_UTILS_H