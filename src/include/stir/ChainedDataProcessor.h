//
//
/*
    Copyright (C) 2000- 2011, Hammersmith Imanet Ltd
    This file is part of STIR.

    SPDX-License-Identifier: Apache-2.0

    See STIR/LICENSE.txt for details
*/
/*!

  \file
  \ingroup DataProcessor
  \brief Declaration of class stir::ChainedDataProcessor

  \author Kris Thielemans

*/
#ifndef __stir_ChainedDataProcessor_H__
#define __stir_ChainedDataProcessor_H__

#include "stir/RegisteredParsingObject.h"
#include "stir/DataProcessor.h"
#include "stir/shared_ptr.h"

START_NAMESPACE_STIR

/*!
  \ingroup DataProcessor
  \brief A class in the DataProcessor hierarchy that calls
   2 DataProcessors in sequence.

  As it is derived from RegisteredParsingObject, it implements all the
  necessary things to parse parameter files etc.
  \warning The 2 argument version of  ChainedDataProcessor::apply
  calls the first data processor with a temporary output data with
  the same characteristics as the input data.


  \warning ChainedDataProcessor::set_up builds only the data
  processor of the first in the data processor chain. This
  is because at this point, we do not really know what the first
  data processor will do to the data (it might change index
  ranges or so), so it is impossible to build the 2nd data
  processor without actually letting the first data processor
  process the data (which might be far too expensive). This is not a real
  problem however, as  ChainedDataProcessor::apply is
  fine as it will  call virtual_set_up for the 2nd
  data processor anyway.

  \par Parameters used for parsing
  \verbatim
  Chained Data Processor Parameters:=
  Data Processor to apply first:=
  Data Processor to apply second:=
  END Chained Data Processor Parameters:=
  \endverbatim
 */

template <typename DataT>
class ChainedDataProcessor
    : public RegisteredParsingObject<ChainedDataProcessor<DataT>, DataProcessor<DataT>, DataProcessor<DataT>>
{
private:
  typedef RegisteredParsingObject<ChainedDataProcessor<DataT>, DataProcessor<DataT>, DataProcessor<DataT>> base_type;

public:
  static const char* const registered_name;

  //! Construct given DataProcessor parameters
  explicit ChainedDataProcessor(shared_ptr<DataProcessor<DataT>> apply_first = shared_ptr<DataProcessor<DataT>>(),
                                shared_ptr<DataProcessor<DataT>> apply_second = shared_ptr<DataProcessor<DataT>>());

private:
  shared_ptr<DataProcessor<DataT>> apply_first;
  shared_ptr<DataProcessor<DataT>> apply_second;

  void set_defaults() override;
  void initialise_keymap() override;

  Succeeded virtual_set_up(const DataT& data) override;

  void virtual_apply(DataT& out_data, const DataT& in_data) const override;
  void virtual_apply(DataT& data) const override;
};

END_NAMESPACE_STIR

#endif
