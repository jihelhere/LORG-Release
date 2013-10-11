// -*- mode: c++ -*-
#pragma once

#include "utils/Tokeniser.h"

class TokeniserFactory {
 public:
  enum TokeniserType {English};

  static Tokeniser * create_tokeniser();
};
