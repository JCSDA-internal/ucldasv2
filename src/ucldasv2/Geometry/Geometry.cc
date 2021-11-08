/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/util/Config.h"

#include "ucldasv2/Geometry/Geometry.h"
#include "ucldasv2/GeometryIterator/GeometryIterator.h"

#include "eckit/config/YAMLConfiguration.h"

#include "oops/util/abor1_cpp.h"

namespace ucldasv2 {

// ----------------------------------------------------------------------------

  Geometry::Geometry(const eckit::Configuration & conf,
                     const eckit::mpi::Comm & comm)
    : comm_(comm),
      fmsinput_(comm, conf) {

    fmsinput_.updateNameList();

    ucldasv2_geo_setup_f90(keyGeom_, &conf, &comm);

    // Set ATLAS lon/lat field
    atlasFieldSet_.reset(new atlas::FieldSet());
    ucldasv2_geo_set_atlas_lonlat_f90(keyGeom_, atlasFieldSet_->get());
    atlas::Field atlasField = atlasFieldSet_->field("lonlat");

    // Create ATLAS function space
    atlasFunctionSpace_.reset(new atlas::functionspace::PointCloud(atlasField));

    // Set ATLAS function space pointer in Fortran
    ucldasv2_geo_set_atlas_functionspace_pointer_f90(keyGeom_,
      atlasFunctionSpace_->get());

    // Fill ATLAS fieldset
    atlasFieldSet_.reset(new atlas::FieldSet());
    ucldasv2_geo_fill_atlas_fieldset_f90(keyGeom_, atlasFieldSet_->get());

  }

// ----------------------------------------------------------------------------

  Geometry::Geometry(const Geometry & other)
    : comm_(other.comm_),
      fmsinput_(other.fmsinput_) {
    const int key_geo = other.keyGeom_;
    ucldasv2_geo_clone_f90(keyGeom_, key_geo);
    atlasFunctionSpace_.reset(new atlas::functionspace::PointCloud(
                              other.atlasFunctionSpace_->lonlat()));
    ucldasv2_geo_set_atlas_functionspace_pointer_f90(keyGeom_,
      atlasFunctionSpace_->get());
    atlasFieldSet_.reset(new atlas::FieldSet());
    for (int jfield = 0; jfield < other.atlasFieldSet_->size(); ++jfield) {
      atlas::Field atlasField = other.atlasFieldSet_->field(jfield);
      atlasFieldSet_->add(atlasField);
    }
  }

// ----------------------------------------------------------------------------

  Geometry::~Geometry() {
    ucldasv2_geo_delete_f90(keyGeom_);
  }

// ----------------------------------------------------------------------------

  GeometryIterator Geometry::begin() const {
    int ist, iend, jst, jend;
    ucldasv2_geo_start_end_f90(keyGeom_, ist, iend, jst, jend);
    return GeometryIterator(*this, ist, jst);
  }

// ----------------------------------------------------------------------------

  GeometryIterator Geometry::end() const {
    // return end of the geometry on this mpi tile
    // decided to return index out of bounds for the iterator loops to work
    return GeometryIterator(*this, -1, -1);
  }

// ----------------------------------------------------------------------------
  void Geometry::print(std::ostream & os) const {
    // TODO(Travis): Implement this correctly.
  }

// ----------------------------------------------------------------------------
  std::vector<double> Geometry::verticalCoord(std::string &) const {
    util::abor1_cpp("Geometry::verticalCoord() needs to be implemented.",
                    __FILE__, __LINE__);
    return {};
  }

  // -----------------------------------------------------------------------------
  atlas::FunctionSpace * Geometry::atlasFunctionSpace() const {
    return atlasFunctionSpace_.get();
  }
  // -----------------------------------------------------------------------------
  atlas::FieldSet * Geometry::atlasFieldSet() const {
    return atlasFieldSet_.get();
  }

// ----------------------------------------------------------------------------

}  // namespace ucldasv2
