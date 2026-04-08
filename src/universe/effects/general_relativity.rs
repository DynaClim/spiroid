use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone)]
pub enum GRModel {
    #[default]
    Disabled,
    Enabled,
}

impl GRModel {
    pub(crate) fn is_enabled(&self) -> bool {
        matches!(self, Self::Enabled)
    }
}
