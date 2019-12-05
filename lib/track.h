#ifndef SYNC_TRACK_H
#define SYNC_TRACK_H

#include <string.h>
#include <stdlib.h>
#include "base.h"

enum key_type {
	KEY_STEP,   /* stay constant */
	KEY_LINEAR, /* lerp to the next value */
	KEY_SMOOTH, /* smooth curve to the next value */
	KEY_RAMP,
	KEY_IN_QUAD,
	KEY_OUT_QUAD,
	KEY_IN_OUT_QUAD,
	KEY_IN_CUBIC,
	KEY_OUT_CUBIC,
	KEY_IN_OUT_CUBIC,
	KEY_IN_QUART,
	KEY_OUT_QUART,
	KEY_IN_OUT_QUART,
	KEY_IN_QUINT,
	KEY_OUT_QUINT,
	KEY_IN_OUT_QUINT,
	KEY_IN_SINE,
	KEY_OUT_SINE,
	KEY_IN_OUT_SINE,
	KEY_IN_CIRC,
	KEY_OUT_CIRC,
	KEY_IN_OUT_CIRC,
	KEY_IN_EXPO,
	KEY_OUT_EXPO,
	KEY_IN_OUT_EXPO,
	KEY_IN_ELASTIC,
	KEY_OUT_ELASTIC,
	KEY_IN_OUT_ELASTIC,
	KEY_IN_BACK,
	KEY_OUT_BACK,
	KEY_IN_OUT_BACK,
	KEY_IN_BOUNCE,
	KEY_OUT_BOUNCE,
	KEY_IN_OUT_BOUNCE,

	KEY_TYPE_COUNT
};

struct track_key {
	int row;
	float value;
	enum key_type type;
};

struct sync_track {
	char *name;
	struct track_key *keys;
	int num_keys;
};

int sync_find_key(const struct sync_track *, int);
static inline int key_idx_floor(const struct sync_track *t, int row)
{
	int idx = sync_find_key(t, row);
	if (idx < 0)
		idx = -idx - 2;
	return idx;
}

#ifndef SYNC_PLAYER
int sync_set_key(struct sync_track *, const struct track_key *);
int sync_del_key(struct sync_track *, int);
static inline int is_key_frame(const struct sync_track *t, int row)
{
	return sync_find_key(t, row) >= 0;
}

#endif /* !defined(SYNC_PLAYER) */

#endif /* SYNC_TRACK_H */
