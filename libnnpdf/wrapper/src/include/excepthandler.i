%exception {
  try {
    $action 
  } 
  catch(NNPDF::FileError &_e) {
    SWIG_exception(SWIG_IOError, const_cast<char*>(_e.what()));
  }
  catch(NNPDF::EvaluationError &_e) {
    SWIG_exception(SWIG_ValueError, const_cast<char*>(_e.what()));
  }
  catch(NNPDF::InitError &_e) {
    SWIG_exception(SWIG_RuntimeError, const_cast<char*>(_e.what()));
  }
  catch(NNPDF::RangeError &_e) {
    SWIG_exception(SWIG_ValueError, const_cast<char*>(_e.what()));
  }
  catch(NNPDF::LengthError &_e) {
    SWIG_exception(SWIG_IndexError, const_cast<char*>(_e.what()));
  }
  catch(NNPDF::LogError &_e) {
    SWIG_exception(SWIG_RuntimeError, const_cast<char*>(_e.what()));
  }
  catch(NNPDF::UserError &_e) {
    SWIG_exception(SWIG_ValueError, const_cast<char*>(_e.what()));
  }
  catch(NNPDF::RuntimeException &_e) {
    SWIG_exception(SWIG_RuntimeError, const_cast<char*>(_e.what()));
  }
  catch(NNPDF::LogicException &_e) {
    SWIG_exception(SWIG_ValueError, const_cast<char*>(_e.what()));
  }
  catch (const std::exception& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
}
