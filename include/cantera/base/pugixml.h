#if CT_USE_SYSTEM_PUGIXML
  #include "pugixml/format.h"
  #include "pugixml/ostream.h"
#else
  #include "cantera/ext/pugixml/pugixml.hpp"
  #include "cantera/ext/pugixml/pugiconfig.hpp"
#endif
