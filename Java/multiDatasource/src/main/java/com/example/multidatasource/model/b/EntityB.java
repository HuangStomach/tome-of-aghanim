
package com.example.multidatasource.model.b;

import lombok.Data;
import javax.persistence.*;

@Data
@Entity
@Table(name = "entity_b")
public class EntityB {
    @Id
    @GeneratedValue(strategy = GenerationType.IDENTITY)
    private Long id;
    
    private String code;
    private Integer value;
}
