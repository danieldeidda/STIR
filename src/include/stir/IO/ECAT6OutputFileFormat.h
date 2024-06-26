//
//
/*
    Copyright (C) 2002-2007, Hammersmith Imanet Ltd
    This file is part of STIR.

    SPDX-License-Identifier: Apache-2.0

    See STIR/LICENSE.txt for details
*/
/*!

  \file
  \ingroup ECAT
  \brief Declaration of class stir::ecat::ecat7::ECAT6OutputFileFormat

  \author Kris Thielemans

*/

#ifndef __stir_IO_ECAT6OutputFileFormat_H__
#define __stir_IO_ECAT6OutputFileFormat_H__

#include "stir/IO/OutputFileFormat.h"
#include "stir/RegisteredParsingObject.h"
// include for namespace macros
#include "stir/IO/stir_ecat_common.h"
#include <string>
#include "stir/deprecated.h"

START_NAMESPACE_STIR
template <int num_dimensions, typename elemT>
class DiscretisedDensity;

START_NAMESPACE_ECAT
START_NAMESPACE_ECAT6

/*!
  \ingroup ECAT
  \brief
  Implementation of OutputFileFormat paradigm for the ECAT6 format.

  \warning Currently output always uses 2-byte signed integers in
  little-endian byte order.

  \deprecated
 */
class STIR_DEPRECATED ECAT6OutputFileFormat : public RegisteredParsingObject<ECAT6OutputFileFormat,
                                                                             OutputFileFormat<DiscretisedDensity<3, float>>,
                                                                             OutputFileFormat<DiscretisedDensity<3, float>>>
{
private:
  typedef RegisteredParsingObject<ECAT6OutputFileFormat,
                                  OutputFileFormat<DiscretisedDensity<3, float>>,
                                  OutputFileFormat<DiscretisedDensity<3, float>>>
      base_type;

public:
  //! Name which will be used when parsing an OutputFileFormat object
  static const char* const registered_name;

  ECAT6OutputFileFormat(const NumericType& = NumericType::SHORT, const ByteOrder& = ByteOrder::native);

  //! Set type of numbers to be used for output
  /*! Currently the return value will always be NumericType::SHORT */
  virtual NumericType set_type_of_numbers(const NumericType&, const bool warn = false);
  //! Set byte order to be used for output
  /*! Currently the return value will always be ByteOrder::LITTLEENDIAN */
  virtual ByteOrder set_byte_order(const ByteOrder&, const bool warn = false);

public:
  std::string default_scanner_name;

protected:
  virtual Succeeded actual_write_to_file(std::string& output_filename, const DiscretisedDensity<3, float>& density) const;
  virtual void set_defaults();
  virtual void initialise_keymap();
  virtual bool post_processing();
};

END_NAMESPACE_ECAT6
END_NAMESPACE_ECAT
END_NAMESPACE_STIR

#endif
