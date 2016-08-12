// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdIsrcdIWCSimRootDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimPmtInfo.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_WCSimRootCherenkovDigiHit(void *p = 0);
   static void *newArray_WCSimRootCherenkovDigiHit(Long_t size, void *p);
   static void delete_WCSimRootCherenkovDigiHit(void *p);
   static void deleteArray_WCSimRootCherenkovDigiHit(void *p);
   static void destruct_WCSimRootCherenkovDigiHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRootCherenkovDigiHit*)
   {
      ::WCSimRootCherenkovDigiHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRootCherenkovDigiHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRootCherenkovDigiHit", ::WCSimRootCherenkovDigiHit::Class_Version(), "WCSimRootEvent.hh", 126,
                  typeid(::WCSimRootCherenkovDigiHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCSimRootCherenkovDigiHit::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRootCherenkovDigiHit) );
      instance.SetNew(&new_WCSimRootCherenkovDigiHit);
      instance.SetNewArray(&newArray_WCSimRootCherenkovDigiHit);
      instance.SetDelete(&delete_WCSimRootCherenkovDigiHit);
      instance.SetDeleteArray(&deleteArray_WCSimRootCherenkovDigiHit);
      instance.SetDestructor(&destruct_WCSimRootCherenkovDigiHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRootCherenkovDigiHit*)
   {
      return GenerateInitInstanceLocal((::WCSimRootCherenkovDigiHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRootCherenkovDigiHit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_WCSimRootCherenkovHit(void *p = 0);
   static void *newArray_WCSimRootCherenkovHit(Long_t size, void *p);
   static void delete_WCSimRootCherenkovHit(void *p);
   static void deleteArray_WCSimRootCherenkovHit(void *p);
   static void destruct_WCSimRootCherenkovHit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRootCherenkovHit*)
   {
      ::WCSimRootCherenkovHit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRootCherenkovHit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRootCherenkovHit", ::WCSimRootCherenkovHit::Class_Version(), "WCSimRootEvent.hh", 84,
                  typeid(::WCSimRootCherenkovHit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCSimRootCherenkovHit::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRootCherenkovHit) );
      instance.SetNew(&new_WCSimRootCherenkovHit);
      instance.SetNewArray(&newArray_WCSimRootCherenkovHit);
      instance.SetDelete(&delete_WCSimRootCherenkovHit);
      instance.SetDeleteArray(&deleteArray_WCSimRootCherenkovHit);
      instance.SetDestructor(&destruct_WCSimRootCherenkovHit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRootCherenkovHit*)
   {
      return GenerateInitInstanceLocal((::WCSimRootCherenkovHit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRootCherenkovHit*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_WCSimRootCherenkovHitTime(void *p = 0);
   static void *newArray_WCSimRootCherenkovHitTime(Long_t size, void *p);
   static void delete_WCSimRootCherenkovHitTime(void *p);
   static void deleteArray_WCSimRootCherenkovHitTime(void *p);
   static void destruct_WCSimRootCherenkovHitTime(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRootCherenkovHitTime*)
   {
      ::WCSimRootCherenkovHitTime *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRootCherenkovHitTime >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRootCherenkovHitTime", ::WCSimRootCherenkovHitTime::Class_Version(), "WCSimRootEvent.hh", 103,
                  typeid(::WCSimRootCherenkovHitTime), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCSimRootCherenkovHitTime::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRootCherenkovHitTime) );
      instance.SetNew(&new_WCSimRootCherenkovHitTime);
      instance.SetNewArray(&newArray_WCSimRootCherenkovHitTime);
      instance.SetDelete(&delete_WCSimRootCherenkovHitTime);
      instance.SetDeleteArray(&deleteArray_WCSimRootCherenkovHitTime);
      instance.SetDestructor(&destruct_WCSimRootCherenkovHitTime);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRootCherenkovHitTime*)
   {
      return GenerateInitInstanceLocal((::WCSimRootCherenkovHitTime*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRootCherenkovHitTime*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_WCSimRootTrack(void *p = 0);
   static void *newArray_WCSimRootTrack(Long_t size, void *p);
   static void delete_WCSimRootTrack(void *p);
   static void deleteArray_WCSimRootTrack(void *p);
   static void destruct_WCSimRootTrack(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRootTrack*)
   {
      ::WCSimRootTrack *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRootTrack >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRootTrack", ::WCSimRootTrack::Class_Version(), "WCSimRootEvent.hh", 24,
                  typeid(::WCSimRootTrack), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCSimRootTrack::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRootTrack) );
      instance.SetNew(&new_WCSimRootTrack);
      instance.SetNewArray(&newArray_WCSimRootTrack);
      instance.SetDelete(&delete_WCSimRootTrack);
      instance.SetDeleteArray(&deleteArray_WCSimRootTrack);
      instance.SetDestructor(&destruct_WCSimRootTrack);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRootTrack*)
   {
      return GenerateInitInstanceLocal((::WCSimRootTrack*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRootTrack*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_WCSimRootEventHeader(void *p = 0);
   static void *newArray_WCSimRootEventHeader(Long_t size, void *p);
   static void delete_WCSimRootEventHeader(void *p);
   static void deleteArray_WCSimRootEventHeader(void *p);
   static void destruct_WCSimRootEventHeader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRootEventHeader*)
   {
      ::WCSimRootEventHeader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRootEventHeader >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRootEventHeader", ::WCSimRootEventHeader::Class_Version(), "WCSimRootEvent.hh", 152,
                  typeid(::WCSimRootEventHeader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCSimRootEventHeader::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRootEventHeader) );
      instance.SetNew(&new_WCSimRootEventHeader);
      instance.SetNewArray(&newArray_WCSimRootEventHeader);
      instance.SetDelete(&delete_WCSimRootEventHeader);
      instance.SetDeleteArray(&deleteArray_WCSimRootEventHeader);
      instance.SetDestructor(&destruct_WCSimRootEventHeader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRootEventHeader*)
   {
      return GenerateInitInstanceLocal((::WCSimRootEventHeader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRootEventHeader*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_WCSimRootTrigger(void *p = 0);
   static void *newArray_WCSimRootTrigger(Long_t size, void *p);
   static void delete_WCSimRootTrigger(void *p);
   static void deleteArray_WCSimRootTrigger(void *p);
   static void destruct_WCSimRootTrigger(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRootTrigger*)
   {
      ::WCSimRootTrigger *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRootTrigger >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRootTrigger", ::WCSimRootTrigger::Class_Version(), "WCSimRootEvent.hh", 205,
                  typeid(::WCSimRootTrigger), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCSimRootTrigger::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRootTrigger) );
      instance.SetNew(&new_WCSimRootTrigger);
      instance.SetNewArray(&newArray_WCSimRootTrigger);
      instance.SetDelete(&delete_WCSimRootTrigger);
      instance.SetDeleteArray(&deleteArray_WCSimRootTrigger);
      instance.SetDestructor(&destruct_WCSimRootTrigger);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRootTrigger*)
   {
      return GenerateInitInstanceLocal((::WCSimRootTrigger*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRootTrigger*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_WCSimRootEvent(void *p = 0);
   static void *newArray_WCSimRootEvent(Long_t size, void *p);
   static void delete_WCSimRootEvent(void *p);
   static void deleteArray_WCSimRootEvent(void *p);
   static void destruct_WCSimRootEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRootEvent*)
   {
      ::WCSimRootEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRootEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRootEvent", ::WCSimRootEvent::Class_Version(), "WCSimRootEvent.hh", 327,
                  typeid(::WCSimRootEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCSimRootEvent::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRootEvent) );
      instance.SetNew(&new_WCSimRootEvent);
      instance.SetNewArray(&newArray_WCSimRootEvent);
      instance.SetDelete(&delete_WCSimRootEvent);
      instance.SetDeleteArray(&deleteArray_WCSimRootEvent);
      instance.SetDestructor(&destruct_WCSimRootEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRootEvent*)
   {
      return GenerateInitInstanceLocal((::WCSimRootEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRootEvent*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_WCSimRootPi0(void *p = 0);
   static void *newArray_WCSimRootPi0(Long_t size, void *p);
   static void delete_WCSimRootPi0(void *p);
   static void deleteArray_WCSimRootPi0(void *p);
   static void destruct_WCSimRootPi0(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRootPi0*)
   {
      ::WCSimRootPi0 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRootPi0 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRootPi0", ::WCSimRootPi0::Class_Version(), "WCSimRootEvent.hh", 176,
                  typeid(::WCSimRootPi0), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCSimRootPi0::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRootPi0) );
      instance.SetNew(&new_WCSimRootPi0);
      instance.SetNewArray(&newArray_WCSimRootPi0);
      instance.SetDelete(&delete_WCSimRootPi0);
      instance.SetDeleteArray(&deleteArray_WCSimRootPi0);
      instance.SetDestructor(&destruct_WCSimRootPi0);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRootPi0*)
   {
      return GenerateInitInstanceLocal((::WCSimRootPi0*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRootPi0*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_WCSimRootGeom(void *p = 0);
   static void *newArray_WCSimRootGeom(Long_t size, void *p);
   static void delete_WCSimRootGeom(void *p);
   static void deleteArray_WCSimRootGeom(void *p);
   static void destruct_WCSimRootGeom(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRootGeom*)
   {
      ::WCSimRootGeom *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRootGeom >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRootGeom", ::WCSimRootGeom::Class_Version(), "WCSimRootGeom.hh", 49,
                  typeid(::WCSimRootGeom), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCSimRootGeom::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRootGeom) );
      instance.SetNew(&new_WCSimRootGeom);
      instance.SetNewArray(&newArray_WCSimRootGeom);
      instance.SetDelete(&delete_WCSimRootGeom);
      instance.SetDeleteArray(&deleteArray_WCSimRootGeom);
      instance.SetDestructor(&destruct_WCSimRootGeom);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRootGeom*)
   {
      return GenerateInitInstanceLocal((::WCSimRootGeom*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRootGeom*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_WCSimRootPMT(void *p = 0);
   static void *newArray_WCSimRootPMT(Long_t size, void *p);
   static void delete_WCSimRootPMT(void *p);
   static void deleteArray_WCSimRootPMT(void *p);
   static void destruct_WCSimRootPMT(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimRootPMT*)
   {
      ::WCSimRootPMT *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimRootPMT >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimRootPMT", ::WCSimRootPMT::Class_Version(), "WCSimRootGeom.hh", 20,
                  typeid(::WCSimRootPMT), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCSimRootPMT::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimRootPMT) );
      instance.SetNew(&new_WCSimRootPMT);
      instance.SetNewArray(&newArray_WCSimRootPMT);
      instance.SetDelete(&delete_WCSimRootPMT);
      instance.SetDeleteArray(&deleteArray_WCSimRootPMT);
      instance.SetDestructor(&destruct_WCSimRootPMT);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimRootPMT*)
   {
      return GenerateInitInstanceLocal((::WCSimRootPMT*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimRootPMT*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_WCSimPmtInfo(void *p = 0);
   static void *newArray_WCSimPmtInfo(Long_t size, void *p);
   static void delete_WCSimPmtInfo(void *p);
   static void deleteArray_WCSimPmtInfo(void *p);
   static void destruct_WCSimPmtInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WCSimPmtInfo*)
   {
      ::WCSimPmtInfo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WCSimPmtInfo >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WCSimPmtInfo", ::WCSimPmtInfo::Class_Version(), "WCSimPmtInfo.hh", 14,
                  typeid(::WCSimPmtInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::WCSimPmtInfo::Dictionary, isa_proxy, 4,
                  sizeof(::WCSimPmtInfo) );
      instance.SetNew(&new_WCSimPmtInfo);
      instance.SetNewArray(&newArray_WCSimPmtInfo);
      instance.SetDelete(&delete_WCSimPmtInfo);
      instance.SetDeleteArray(&deleteArray_WCSimPmtInfo);
      instance.SetDestructor(&destruct_WCSimPmtInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WCSimPmtInfo*)
   {
      return GenerateInitInstanceLocal((::WCSimPmtInfo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WCSimPmtInfo*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr WCSimRootCherenkovDigiHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRootCherenkovDigiHit::Class_Name()
{
   return "WCSimRootCherenkovDigiHit";
}

//______________________________________________________________________________
const char *WCSimRootCherenkovDigiHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootCherenkovDigiHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRootCherenkovDigiHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootCherenkovDigiHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCSimRootCherenkovDigiHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootCherenkovDigiHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRootCherenkovDigiHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootCherenkovDigiHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WCSimRootCherenkovHit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRootCherenkovHit::Class_Name()
{
   return "WCSimRootCherenkovHit";
}

//______________________________________________________________________________
const char *WCSimRootCherenkovHit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootCherenkovHit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRootCherenkovHit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootCherenkovHit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCSimRootCherenkovHit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootCherenkovHit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRootCherenkovHit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootCherenkovHit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WCSimRootCherenkovHitTime::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRootCherenkovHitTime::Class_Name()
{
   return "WCSimRootCherenkovHitTime";
}

//______________________________________________________________________________
const char *WCSimRootCherenkovHitTime::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootCherenkovHitTime*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRootCherenkovHitTime::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootCherenkovHitTime*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCSimRootCherenkovHitTime::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootCherenkovHitTime*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRootCherenkovHitTime::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootCherenkovHitTime*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WCSimRootTrack::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRootTrack::Class_Name()
{
   return "WCSimRootTrack";
}

//______________________________________________________________________________
const char *WCSimRootTrack::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootTrack*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRootTrack::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootTrack*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCSimRootTrack::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootTrack*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRootTrack::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootTrack*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WCSimRootEventHeader::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRootEventHeader::Class_Name()
{
   return "WCSimRootEventHeader";
}

//______________________________________________________________________________
const char *WCSimRootEventHeader::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootEventHeader*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRootEventHeader::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootEventHeader*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCSimRootEventHeader::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootEventHeader*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRootEventHeader::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootEventHeader*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WCSimRootTrigger::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRootTrigger::Class_Name()
{
   return "WCSimRootTrigger";
}

//______________________________________________________________________________
const char *WCSimRootTrigger::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootTrigger*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRootTrigger::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootTrigger*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCSimRootTrigger::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootTrigger*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRootTrigger::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootTrigger*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WCSimRootEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRootEvent::Class_Name()
{
   return "WCSimRootEvent";
}

//______________________________________________________________________________
const char *WCSimRootEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRootEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCSimRootEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRootEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WCSimRootPi0::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRootPi0::Class_Name()
{
   return "WCSimRootPi0";
}

//______________________________________________________________________________
const char *WCSimRootPi0::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootPi0*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRootPi0::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootPi0*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCSimRootPi0::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootPi0*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRootPi0::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootPi0*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WCSimRootGeom::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRootGeom::Class_Name()
{
   return "WCSimRootGeom";
}

//______________________________________________________________________________
const char *WCSimRootGeom::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootGeom*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRootGeom::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootGeom*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCSimRootGeom::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootGeom*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRootGeom::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootGeom*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WCSimRootPMT::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimRootPMT::Class_Name()
{
   return "WCSimRootPMT";
}

//______________________________________________________________________________
const char *WCSimRootPMT::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootPMT*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimRootPMT::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootPMT*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCSimRootPMT::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootPMT*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimRootPMT::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimRootPMT*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WCSimPmtInfo::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WCSimPmtInfo::Class_Name()
{
   return "WCSimPmtInfo";
}

//______________________________________________________________________________
const char *WCSimPmtInfo::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimPmtInfo*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WCSimPmtInfo::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WCSimPmtInfo*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WCSimPmtInfo::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimPmtInfo*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WCSimPmtInfo::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WCSimPmtInfo*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void WCSimRootCherenkovDigiHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRootCherenkovDigiHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRootCherenkovDigiHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRootCherenkovDigiHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimRootCherenkovDigiHit(void *p) {
      return  p ? new(p) ::WCSimRootCherenkovDigiHit : new ::WCSimRootCherenkovDigiHit;
   }
   static void *newArray_WCSimRootCherenkovDigiHit(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimRootCherenkovDigiHit[nElements] : new ::WCSimRootCherenkovDigiHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimRootCherenkovDigiHit(void *p) {
      delete ((::WCSimRootCherenkovDigiHit*)p);
   }
   static void deleteArray_WCSimRootCherenkovDigiHit(void *p) {
      delete [] ((::WCSimRootCherenkovDigiHit*)p);
   }
   static void destruct_WCSimRootCherenkovDigiHit(void *p) {
      typedef ::WCSimRootCherenkovDigiHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRootCherenkovDigiHit

//______________________________________________________________________________
void WCSimRootCherenkovHit::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRootCherenkovHit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRootCherenkovHit::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRootCherenkovHit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimRootCherenkovHit(void *p) {
      return  p ? new(p) ::WCSimRootCherenkovHit : new ::WCSimRootCherenkovHit;
   }
   static void *newArray_WCSimRootCherenkovHit(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimRootCherenkovHit[nElements] : new ::WCSimRootCherenkovHit[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimRootCherenkovHit(void *p) {
      delete ((::WCSimRootCherenkovHit*)p);
   }
   static void deleteArray_WCSimRootCherenkovHit(void *p) {
      delete [] ((::WCSimRootCherenkovHit*)p);
   }
   static void destruct_WCSimRootCherenkovHit(void *p) {
      typedef ::WCSimRootCherenkovHit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRootCherenkovHit

//______________________________________________________________________________
void WCSimRootCherenkovHitTime::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRootCherenkovHitTime.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRootCherenkovHitTime::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRootCherenkovHitTime::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimRootCherenkovHitTime(void *p) {
      return  p ? new(p) ::WCSimRootCherenkovHitTime : new ::WCSimRootCherenkovHitTime;
   }
   static void *newArray_WCSimRootCherenkovHitTime(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimRootCherenkovHitTime[nElements] : new ::WCSimRootCherenkovHitTime[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimRootCherenkovHitTime(void *p) {
      delete ((::WCSimRootCherenkovHitTime*)p);
   }
   static void deleteArray_WCSimRootCherenkovHitTime(void *p) {
      delete [] ((::WCSimRootCherenkovHitTime*)p);
   }
   static void destruct_WCSimRootCherenkovHitTime(void *p) {
      typedef ::WCSimRootCherenkovHitTime current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRootCherenkovHitTime

//______________________________________________________________________________
void WCSimRootTrack::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRootTrack.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRootTrack::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRootTrack::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimRootTrack(void *p) {
      return  p ? new(p) ::WCSimRootTrack : new ::WCSimRootTrack;
   }
   static void *newArray_WCSimRootTrack(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimRootTrack[nElements] : new ::WCSimRootTrack[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimRootTrack(void *p) {
      delete ((::WCSimRootTrack*)p);
   }
   static void deleteArray_WCSimRootTrack(void *p) {
      delete [] ((::WCSimRootTrack*)p);
   }
   static void destruct_WCSimRootTrack(void *p) {
      typedef ::WCSimRootTrack current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRootTrack

//______________________________________________________________________________
void WCSimRootEventHeader::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRootEventHeader.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRootEventHeader::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRootEventHeader::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimRootEventHeader(void *p) {
      return  p ? new(p) ::WCSimRootEventHeader : new ::WCSimRootEventHeader;
   }
   static void *newArray_WCSimRootEventHeader(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimRootEventHeader[nElements] : new ::WCSimRootEventHeader[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimRootEventHeader(void *p) {
      delete ((::WCSimRootEventHeader*)p);
   }
   static void deleteArray_WCSimRootEventHeader(void *p) {
      delete [] ((::WCSimRootEventHeader*)p);
   }
   static void destruct_WCSimRootEventHeader(void *p) {
      typedef ::WCSimRootEventHeader current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRootEventHeader

//______________________________________________________________________________
void WCSimRootTrigger::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRootTrigger.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRootTrigger::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRootTrigger::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimRootTrigger(void *p) {
      return  p ? new(p) ::WCSimRootTrigger : new ::WCSimRootTrigger;
   }
   static void *newArray_WCSimRootTrigger(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimRootTrigger[nElements] : new ::WCSimRootTrigger[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimRootTrigger(void *p) {
      delete ((::WCSimRootTrigger*)p);
   }
   static void deleteArray_WCSimRootTrigger(void *p) {
      delete [] ((::WCSimRootTrigger*)p);
   }
   static void destruct_WCSimRootTrigger(void *p) {
      typedef ::WCSimRootTrigger current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRootTrigger

//______________________________________________________________________________
void WCSimRootEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRootEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRootEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRootEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimRootEvent(void *p) {
      return  p ? new(p) ::WCSimRootEvent : new ::WCSimRootEvent;
   }
   static void *newArray_WCSimRootEvent(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimRootEvent[nElements] : new ::WCSimRootEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimRootEvent(void *p) {
      delete ((::WCSimRootEvent*)p);
   }
   static void deleteArray_WCSimRootEvent(void *p) {
      delete [] ((::WCSimRootEvent*)p);
   }
   static void destruct_WCSimRootEvent(void *p) {
      typedef ::WCSimRootEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRootEvent

//______________________________________________________________________________
void WCSimRootPi0::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRootPi0.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRootPi0::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRootPi0::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimRootPi0(void *p) {
      return  p ? new(p) ::WCSimRootPi0 : new ::WCSimRootPi0;
   }
   static void *newArray_WCSimRootPi0(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimRootPi0[nElements] : new ::WCSimRootPi0[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimRootPi0(void *p) {
      delete ((::WCSimRootPi0*)p);
   }
   static void deleteArray_WCSimRootPi0(void *p) {
      delete [] ((::WCSimRootPi0*)p);
   }
   static void destruct_WCSimRootPi0(void *p) {
      typedef ::WCSimRootPi0 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRootPi0

//______________________________________________________________________________
void WCSimRootGeom::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRootGeom.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRootGeom::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRootGeom::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimRootGeom(void *p) {
      return  p ? new(p) ::WCSimRootGeom : new ::WCSimRootGeom;
   }
   static void *newArray_WCSimRootGeom(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimRootGeom[nElements] : new ::WCSimRootGeom[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimRootGeom(void *p) {
      delete ((::WCSimRootGeom*)p);
   }
   static void deleteArray_WCSimRootGeom(void *p) {
      delete [] ((::WCSimRootGeom*)p);
   }
   static void destruct_WCSimRootGeom(void *p) {
      typedef ::WCSimRootGeom current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRootGeom

//______________________________________________________________________________
void WCSimRootPMT::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimRootPMT.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimRootPMT::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimRootPMT::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimRootPMT(void *p) {
      return  p ? new(p) ::WCSimRootPMT : new ::WCSimRootPMT;
   }
   static void *newArray_WCSimRootPMT(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimRootPMT[nElements] : new ::WCSimRootPMT[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimRootPMT(void *p) {
      delete ((::WCSimRootPMT*)p);
   }
   static void deleteArray_WCSimRootPMT(void *p) {
      delete [] ((::WCSimRootPMT*)p);
   }
   static void destruct_WCSimRootPMT(void *p) {
      typedef ::WCSimRootPMT current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimRootPMT

//______________________________________________________________________________
void WCSimPmtInfo::Streamer(TBuffer &R__b)
{
   // Stream an object of class WCSimPmtInfo.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WCSimPmtInfo::Class(),this);
   } else {
      R__b.WriteClassBuffer(WCSimPmtInfo::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_WCSimPmtInfo(void *p) {
      return  p ? new(p) ::WCSimPmtInfo : new ::WCSimPmtInfo;
   }
   static void *newArray_WCSimPmtInfo(Long_t nElements, void *p) {
      return p ? new(p) ::WCSimPmtInfo[nElements] : new ::WCSimPmtInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_WCSimPmtInfo(void *p) {
      delete ((::WCSimPmtInfo*)p);
   }
   static void deleteArray_WCSimPmtInfo(void *p) {
      delete [] ((::WCSimPmtInfo*)p);
   }
   static void destruct_WCSimPmtInfo(void *p) {
      typedef ::WCSimPmtInfo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WCSimPmtInfo

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 214,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEfloatgR_Dictionary();
   static void vectorlEfloatgR_TClassManip(TClass*);
   static void *new_vectorlEfloatgR(void *p = 0);
   static void *newArray_vectorlEfloatgR(Long_t size, void *p);
   static void delete_vectorlEfloatgR(void *p);
   static void deleteArray_vectorlEfloatgR(void *p);
   static void destruct_vectorlEfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<float>*)
   {
      vector<float> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<float>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<float>", -2, "vector", 214,
                  typeid(vector<float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<float>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<float>*)0x0)->GetClass();
      vectorlEfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEfloatgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float> : new vector<float>;
   }
   static void *newArray_vectorlEfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float>[nElements] : new vector<float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEfloatgR(void *p) {
      delete ((vector<float>*)p);
   }
   static void deleteArray_vectorlEfloatgR(void *p) {
      delete [] ((vector<float>*)p);
   }
   static void destruct_vectorlEfloatgR(void *p) {
      typedef vector<float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<float>

namespace {
  void TriggerDictionaryInitialization_WCSimRootDict_Impl() {
    static const char* headers[] = {
"WCSimRootEvent.hh",
"WCSimRootGeom.hh",
"WCSimPmtInfo.hh",
0
    };
    static const char* includePaths[] = {
"./include",
"/home/marc/LinuxSystemFiles/ROOT/root-6.06.04/build/include",
"/home/marc/LinuxSystemFiles/ROOT/root-6.06.04/build/include",
"/home/marc/LinuxSystemFiles/WCSim/WCSim_v0/gitver/WCSim/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "WCSimRootDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$WCSimRootEvent.hh")))  WCSimRootCherenkovDigiHit;
class __attribute__((annotate("$clingAutoload$WCSimRootEvent.hh")))  WCSimRootCherenkovHit;
class __attribute__((annotate("$clingAutoload$WCSimRootEvent.hh")))  WCSimRootCherenkovHitTime;
class __attribute__((annotate("$clingAutoload$WCSimRootEvent.hh")))  WCSimRootTrack;
class __attribute__((annotate(R"ATTRDUMP(WCSimRootEvent Header)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$WCSimRootEvent.hh")))  WCSimRootEventHeader;
class __attribute__((annotate(R"ATTRDUMP(WCSimRootEvent structure)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$WCSimRootEvent.hh")))  WCSimRootTrigger;
class __attribute__((annotate("$clingAutoload$WCSimRootEvent.hh")))  WCSimRootEvent;
class __attribute__((annotate("$clingAutoload$WCSimRootEvent.hh")))  WCSimRootPi0;
class __attribute__((annotate(R"ATTRDUMP(WCSimRootEvent structure)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$WCSimRootGeom.hh")))  WCSimRootGeom;
class __attribute__((annotate(R"ATTRDUMP(WCSimPMT structure)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$WCSimRootGeom.hh")))  WCSimRootPMT;
class __attribute__((annotate("$clingAutoload$WCSimPmtInfo.hh")))  WCSimPmtInfo;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "WCSimRootDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimPmtInfo.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"WCSimPmtInfo", payloadCode, "@",
"WCSimRootCherenkovDigiHit", payloadCode, "@",
"WCSimRootCherenkovHit", payloadCode, "@",
"WCSimRootCherenkovHitTime", payloadCode, "@",
"WCSimRootEvent", payloadCode, "@",
"WCSimRootEventHeader", payloadCode, "@",
"WCSimRootGeom", payloadCode, "@",
"WCSimRootPMT", payloadCode, "@",
"WCSimRootPi0", payloadCode, "@",
"WCSimRootTrack", payloadCode, "@",
"WCSimRootTrigger", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("WCSimRootDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_WCSimRootDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_WCSimRootDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_WCSimRootDict() {
  TriggerDictionaryInitialization_WCSimRootDict_Impl();
}
