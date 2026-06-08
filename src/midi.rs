use crate::wave::WaveState;

/// A single MIDI event derived from a wave pattern.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct MidiEvent {
    /// MIDI note number (0–127). Mapped from node position.
    pub note: u8,
    /// Velocity (0–127). Mapped from wave amplitude.
    pub velocity: u8,
    /// Time offset in ticks.
    pub tick: u32,
    /// Duration in ticks.
    pub duration: u32,
}

/// Map a wave state to a sequence of MIDI events.
///
/// - **Position → pitch**: node index is linearly mapped to MIDI note range.
/// - **Amplitude → velocity**: displacement magnitude is mapped to velocity.
/// - **Frequency → tempo**: optionally controls the tick spacing.
pub fn wave_to_midi(
    state: &WaveState,
    min_note: u8,
    max_note: u8,
    ticks_per_event: u32,
    duration_ticks: u32,
) -> Vec<MidiEvent> {
    let n = state.displacement.len();
    if n == 0 || min_note >= max_note {
        return vec![];
    }

    let max_amp = state.displacement.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
    if max_amp < 1e-12 {
        return vec![];
    }

    let note_range = (max_note - min_note) as f64;

    state
        .displacement
        .iter()
        .enumerate()
        .filter(|(_, &amp)| amp.abs() > max_amp * 0.05)
        .map(|(i, &amp)| {
            let note_frac = if n > 1 { i as f64 / (n - 1) as f64 } else { 0.5 };
            let note = min_note + (note_frac * note_range).round() as u8;
            let velocity = ((amp.abs() / max_amp) * 127.0).round() as u8;
            MidiEvent {
                note,
                velocity: velocity.clamp(1, 127),
                tick: state.time.round() as u32 + i as u32 * ticks_per_event,
                duration: duration_ticks,
            }
        })
        .collect()
}

/// Convert MIDI events to a simple text-based MIDI file representation (for debugging).
#[allow(dead_code)]
pub fn midi_to_text(events: &[MidiEvent]) -> String {
    events
        .iter()
        .map(|e| format!("tick={:04} note={:03} vel={:03} dur={}", e.tick, e.note, e.velocity, e.duration))
        .collect::<Vec<_>>()
        .join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wave_to_midi_basic() {
        let state = WaveState {
            displacement: vec![0.0, 0.5, 1.0, 0.5, 0.0],
            velocity: vec![0.0; 5],
            time: 0.0,
        };
        let events = wave_to_midi(&state, 60, 72, 120, 480);
        // Nodes 0 and 4 are near zero, should be filtered
        assert!(events.len() >= 2);
        assert!(events.len() <= 5);
        // All notes in range
        for e in &events {
            assert!(e.note >= 60 && e.note <= 72);
            assert!(e.velocity > 0);
        }
    }

    #[test]
    fn test_midi_velocity_scales() {
        let state = WaveState {
            displacement: vec![0.1, 0.5, 1.0],
            velocity: vec![0.0; 3],
            time: 0.0,
        };
        let events = wave_to_midi(&state, 60, 72, 120, 480);
        // Node 2 (amp 1.0) should have highest velocity
        let max_vel_event = events.iter().max_by_key(|e| e.velocity).unwrap();
        assert_eq!(max_vel_event.note, 72); // highest position = max_note
    }

    #[test]
    fn test_empty_wave() {
        let state = WaveState {
            displacement: vec![],
            velocity: vec![],
            time: 0.0,
        };
        let events = wave_to_midi(&state, 60, 72, 120, 480);
        assert!(events.is_empty());
    }

    #[test]
    fn test_zero_amplitude() {
        let state = WaveState {
            displacement: vec![0.0, 0.0, 0.0],
            velocity: vec![0.0; 3],
            time: 0.0,
        };
        let events = wave_to_midi(&state, 60, 72, 120, 480);
        assert!(events.is_empty());
    }

    #[test]
    fn test_midi_to_text() {
        let events = vec![MidiEvent {
            note: 60,
            velocity: 100,
            tick: 0,
            duration: 480,
        }];
        let text = midi_to_text(&events);
        assert!(text.contains("note=060"));
        assert!(text.contains("vel=100"));
    }

    #[test]
    fn test_midi_note_range_single_node() {
        let state = WaveState {
            displacement: vec![1.0],
            velocity: vec![0.0],
            time: 0.0,
        };
        let events = wave_to_midi(&state, 60, 72, 120, 480);
        assert_eq!(events.len(), 1);
        // Single node maps to midpoint
        assert!(events[0].note >= 60);
    }
}
