#ifndef __stir_IO_ITKInputFileFormat_h__
#define __stir_IO_ITKInputFileFormat_h__
/*
    Copyright (C) 2013, Institute for Bioengineering of Catalonia
    Copyright (C) 2014, University College London
    Copyright (C) 2018, Commonwealth Scientific and Industrial Research Organisation
                        Australian eHealth Research Centre
    This file is part of STIR.
    SPDX-License-Identifier: Apache-2.0

    See STIR/LICENSE.txt for details
*/
/*!
  \file
  \ingroup IO
  \brief Declaration of class stir::ITKInputFileFormat

  \author Berta Marti Fuster
  \author Kris Thielemans
  \author Ashley Gillman
*/
#include "stir/IO/InputFileFormat.h"
#include "stir/DiscretisedDensity.h"

START_NAMESPACE_STIR

//! Class for reading images using ITK.
/*! \ingroup IO

    ITK (http://www.itk.org) has its own registry of file formats, so the current class
    provides an interface to that code. We use ITK for reading, and then translate the ITK
    data and meta-info to STIR.

    ITK can read many file formats, see http://www.itk.org/Wiki/ITK/File_Formats for some info.

    This STIR class has special handling for DICOM images. For many modalities, DICOM stores
    each slice in a different file. Normally, ITK reads only a single DICOM file, and hence a single slice.
    As this is not useful for STIR, we use \c itk::GDCMSeriesFileNames to find
    other slices belonging to the same series/time frame/gate as the input filename to read_from_file().

    \warning This translation currently ignores orientation and direction (e.g. of slice order).
*/
template <typename STIRImageType = DiscretisedDensity<3, float>>
class ITKImageInputFileFormat : public InputFileFormat<STIRImageType>
{

public:
  //! This function always returns \c false as ITK cannot read from istream
  bool can_read(const FileSignature& signature, std::istream& input) const override;
  //! Use ITK reader to check if it can read the file (by seeing if it throws an exception or not)
  bool can_read(const FileSignature& signature, const std::string& filename) const override;

  const std::string get_name() const override { return "ITK"; }

protected:
  bool actual_can_read(const FileSignature& signature, std::istream& input) const override;

  //! This function always calls error() as ITK cannot read from istream
  unique_ptr<STIRImageType> read_from_file(std::istream& input) const override;

  //! This function uses ITK for reading and does the translation to STIR
  unique_ptr<STIRImageType> read_from_file(const std::string& filename) const override;
};
END_NAMESPACE_STIR

#endif
