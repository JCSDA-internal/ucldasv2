/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <iomanip>
#include <vector>

#include "ucldasv2/Geometry/Geometry.h"
#include "ucldasv2/Increment/Increment.h"
#include "ucldasv2/State/State.h"
#include "ucldasv2/State/StateFortran.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "ufo/GeoVaLs.h"
#include "ufo/Locations.h"

using oops::Log;

namespace ucldasv2 {

  // ----------------------------------------------------------------------------
  /// Constructor, destructor
  // ----------------------------------------------------------------------------
  State::State(const Geometry & geom, const oops::Variables & vars,
               const util::DateTime & vt)
    : time_(vt), vars_(vars), geom_(new Geometry(geom))
  {
    ucldasv2_state_create_f90(keyFlds_, geom_->toFortran(), vars_);
    Log::trace() << "State::State created." << std::endl;
  }
  // ----------------------------------------------------------------------------
  State::State(const Geometry & geom, const eckit::Configuration & conf)
    : time_(),
      vars_(conf, "state variables"),
      geom_(new Geometry(geom))
  {
    util::DateTime * dtp = &time_;
    oops::Variables vars(vars_);
    ucldasv2_state_create_f90(keyFlds_, geom_->toFortran(), vars);

    if (conf.has("analytic init")) {
      std::string dt;
      conf.get("date", dt);
      time_ = util::DateTime(dt);
      //ucldasv2_state_analytic_f90(toFortran(), &conf, &dtp);
    } else {
      ucldasv2_state_read_file_f90(toFortran(), &conf, &dtp);
    }
    Log::trace() << "State::State created and read in." << std::endl;

  }
  // ----------------------------------------------------------------------------
  State::State(const Geometry & geom, const State & other)
    : vars_(other.vars_), time_(other.time_), geom_(new Geometry(geom))
  {
    ucldasv2_state_create_f90(keyFlds_, geom_->toFortran(), vars_);
    ucldasv2_state_change_resol_f90(toFortran(), other.keyFlds_);
    Log::trace() << "State::State created by interpolation." << std::endl;
  }
  // ----------------------------------------------------------------------------
  State::State(const State & other)
    : vars_(other.vars_), time_(other.time_), geom_(new Geometry(*other.geom_))
  {
    ucldasv2_state_create_f90(keyFlds_, geom_->toFortran(), vars_);
    ucldasv2_state_copy_f90(toFortran(), other.toFortran());
    Log::trace() << "State::State copied." << std::endl;
  }
  // ----------------------------------------------------------------------------
  State::~State() {
    ucldasv2_state_delete_f90(toFortran());
    Log::trace() << "State::State destructed." << std::endl;
  }
  // ----------------------------------------------------------------------------
  State & State::operator+=(const Increment & dx)
  {
    ASSERT(validTime() == dx.validTime());
    // Interpolate increment to analysis grid
    Increment dx_hr(*geom_, dx);

    // Add increment to background state
    ucldasv2_state_add_incr_f90(toFortran(), dx_hr.toFortran());
    return *this;
  }
  // ----------------------------------------------------------------------------
  void State::accumul(const double & zz, const State & xx) {
    ucldasv2_state_axpy_f90(toFortran(), zz, xx.toFortran());
  }
  // ----------------------------------------------------------------------------

  double State::norm() const {
    double zz = 0.0;
    ucldasv2_state_rms_f90(toFortran(), zz);
    return zz;
  }

// ----------------------------------------------------------------------------

  void State::zero() {
    ucldasv2_state_zero_f90(toFortran());
  }

// ----------------------------------------------------------------------------

  void State::read(const eckit::Configuration & files) {
    Log::trace() << "State::State read started." << std::endl;
    util::DateTime * dtp = &time_;
    ucldasv2_state_read_file_f90(toFortran(), &files, &dtp);
    Log::trace() << "State::State read done." << std::endl;
  }

// ----------------------------------------------------------------------------

  void State::write(const eckit::Configuration & files) const {
    const util::DateTime * dtp = &time_;
    ucldasv2_state_write_file_f90(toFortran(), &files, &dtp);
  }

// ----------------------------------------------------------------------------

  size_t State::serialSize() const {
    // Field
    size_t nn;
    ucldasv2_state_serial_size_f90(toFortran(), geom_->toFortran(), nn);

    // Magic factor
    nn += 1;

    // Date and time
    nn += time_.serialSize();
    return nn;
  }
  // -----------------------------------------------------------------------------
  constexpr double SerializeCheckValue = -54321.98765;
  void State::serialize(std::vector<double> & vect) const {
    // Serialize the field
    size_t nn;
    ucldasv2_state_serial_size_f90(toFortran(), geom_->toFortran(), nn);
    std::vector<double> vect_field(nn, 0);
    vect.reserve(vect.size() + nn + 1 + time_.serialSize());
    ucldasv2_state_serialize_f90(toFortran(), geom_->toFortran(), nn,
                             vect_field.data());
    vect.insert(vect.end(), vect_field.begin(), vect_field.end());

    // Magic value placed in serialization; used to validate deserialization
    vect.push_back(SerializeCheckValue);

    // Serialize the date and time
    time_.serialize(vect);
  }
  // ----------------------------------------------------------------------------
  void State::deserialize(const std::vector<double> & vect, size_t & index) {
    // Deserialize the field
    ucldasv2_state_deserialize_f90(toFortran(), geom_->toFortran(), vect.size(),
                               vect.data(), index);

    // Use magic value to validate deserialization
    ASSERT(vect.at(index) == SerializeCheckValue);
    ++index;

    // Deserialize the date and time
    time_.deserialize(vect, index);
  }
  // ----------------------------------------------------------------------------

  std::shared_ptr<const Geometry> State::geometry() const {return geom_;}

  // ----------------------------------------------------------------------------

  void State::print(std::ostream & os) const {
    os << std::endl << "  Valid time: " << validTime();
    int n0, nf;
    ucldasv2_state_sizes_f90(toFortran(), n0, n0, n0, nf);
    std::vector<double> zstat(3*nf);
    ucldasv2_state_gpnorm_f90(toFortran(), nf, zstat[0]);
    for (int jj = 0; jj < nf; ++jj) {
      os << std::endl << std::right << std::setw(7) << vars_[jj]
         << "   min="  <<  std::fixed << std::setw(12) <<
                           std::right << zstat[3*jj]
         << "   max="  <<  std::fixed << std::setw(12) <<
                           std::right << zstat[3*jj+1]
         << "   mean=" <<  std::fixed << std::setw(12) <<
                           std::right << zstat[3*jj+2];
    }
  }

// ----------------------------------------------------------------------------

}  // namespace ucldasv2
