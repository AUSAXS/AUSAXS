#pragma once

class Decomposition {
  public: 
    Decomposition() {}
    virtual ~Decomposition() {}
    virtual void decompose() = 0;
    virtual double determinant() const = 0;

  protected:
    static constexpr double precision = 1e-9;
};