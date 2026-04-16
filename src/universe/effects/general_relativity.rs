use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, PartialEq, Debug, Default, Clone)]
pub enum GeneralRelativityModel {
    #[default]
    Disabled,
    Enabled,
}

impl GeneralRelativityModel {
    pub(crate) fn is_enabled(&self) -> bool {
        matches!(self, Self::Enabled)
    }
}
