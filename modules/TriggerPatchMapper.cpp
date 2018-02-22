#include "TriggerPatchMapper.h"

shared_ptr<TriggerPatchMapper> TriggerPatchMapper::instance = nullptr;

int TriggerPatchMapper::findTriggerPatch( int _channel ){
	if ( nullptr == TriggerPatchMapper::instance ){
		TriggerPatchMapper::instance = shared_ptr<TriggerPatchMapper>( new TriggerPatchMapper() );
	}

	return TriggerPatchMapper::instance->triggerPatch( _channel );
}

int TriggerPatchMapper::findTriggerPatch( int _bl, int _mod ){
	if ( nullptr == TriggerPatchMapper::instance ){
		TriggerPatchMapper::instance = shared_ptr<TriggerPatchMapper>( new TriggerPatchMapper() );
	}

	return TriggerPatchMapper::instance->triggerPatch( _bl, _mod );
}